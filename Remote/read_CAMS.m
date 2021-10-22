function MC = read_CAMS(file,varargin)
% S = READ_CAMS(FILE) - Read a CAMS/McClear data FILE, as returned by WPS_CAMS_MCCLEAR or 
%   manually downloaded from <http://www.soda-pro.com/web-services/radiation/...>
% 
%     # Coding: utf-8
%     # File format version: 3
%     # Title: CAMS McClear v3.1 model of clear-sky irradiation / CAMS Radiation Service v3.2 all-sky irradiation
%     # ...
%     # Date begin (ISO 8601): 2018-01-01T00:00:00.0
%     # Date end (ISO 8601): 2019-01-01T00:00:00.0
%     # Latitude (positive North, ISO 19115): 53.1400
%     #...
%     # Columns:
%     # 1. Observation period (ISO 8601)
%     # 2. TOA. Irradiation on horizontal plane at the top of atmosphere (Wh/m²)
%     # 3. Clear sky GHI. Clear sky global irradiation on horizontal plane at ground level (Wh/m²)
%     # ...
%
% INPUT: FILE name, parsed by PICKFILE.
% OUTPUT: Note change of units in some variables (notably Wh/m² to W/m² in irradiance)
%
%     S.info: cellstring, row headers not parsed already to other fields, without # markers
%     S.location: structure with fields {latitude,longitude,altitude, and CAMS_altitude}
%           (degrees north, degrees east, and mASL, respectively).
%     S.date_begin: [1,6] datevec
%     S.date_end: [1,6] datevec
%     S.time_ref: 'Universal time (UT)'
%     S.dt: (DURATION scalar) time-step
%
%     S.t: [N×1] datenum vector of UTC time-steps (center of interval)
%     S.TOA. Irradiation on horizontal plane at the top of atmosphere (W/m²)
%     S.CSGHI. Clear sky global irradiation on horizontal plane at ground level (W/m²)
%     S.CSBHI. Clear sky beam irradiation on horizontal plane at ground level (W/m²)
%     S.CSDHI. Clear sky diffuse irradiation on horizontal plane at ground level (W/m²)
%     S.CSBNI. Clear sky beam irradiation on mobile plane following the sun at normal incidence (W/m²)
%     S.sza. Solar zenithal angle for the middle of the interval (deg)
%     S.tpw. Total column content of water vapour (kg/m² ~ mm)
%     S.fiso/fvol/fgeo. MODIS-like BRDF parameters fiso/fvol/fgeo
%     S.albedo. Ground albedo
%
% S = READ_CAMS(FILE,VARARGIN) - pass additional arguments to READTXTFILE, e.g....
% S = READ_CAMS(FILE,'data',false) - reads and returns only header-data information.
%
% SEE ALSO: WPS_MERRA2, WPS_CAMS_MCCLEAR, GETMETEODATA, COMPLETEMETEODATA

    if nargin < 1, file = '*'; end
    file = pickfile(file);

    % Updated for version 3.1, 5.02.2019
    VAR_LIST_IN = {'GHI','BNI','DHI','TOA','ClearSkyGHI','ClearSkyBHI','ClearSkyDHI','ClearSkyBNI','sza','fiso','fvol','fgeo','albedo','tcwv'};
    VAR_LIST_OUT = {'GHI','BNI','DHI','TOA','CSGHI','CSBHI','CSDHI','CSBNI','sza','fiso','fvol','fgeo','albedo','tpw'};   
    UNIT_SCALE = @(dt) [60/dt,60/dt,60/dt,60/dt,60/dt,60/dt,60/dt,60/dt,1,1,1,1,1,1];
    
    OPTIONAL = false(size(VAR_LIST_IN));
    OPTIONAL(1:3) = true;       % {GHI,BNI,DHI} missing for McClear files
    OPTIONAL(end-5:end) = true; % {sza, fiso .. albedo} missing for CAMS v2 files
        
    TZ = ''; % should be set by MC.time_ref before first call to @TZdate (§)
    
    % # PARAMNAME;VALUE pairs in file header
    PARAMS = {'Latitude.*','latitude',@str2double;
              'Longitude.*','longitude',@str2double;
              'Altitude.*','altitude',@str2double;
              'Elevation of CAMS cell.*','CAMS_altitude',@str2double; % version >= 3
              'Time reference','time_ref',@checktimezone; % sets persistent TZ (§)...
              'Date begin.*','date_begin',@TZdate; % (§)
              'Date end.*','date_end',@TZdate; %(§)
              'Summarization.*','dt',@parsetimestep};
          
    checkfile(file); % Check if file is not there, empty, or just an XML error code

    data = readtxtfile(file,PARAMS(:,1:2),'regexpmatch',true,'treatasempty','nan',varargin{:});
    MC.info = data.comments;
    
    for j = 1:size(PARAMS,1)
        f = PARAMS{j,2};
        if ~isfield(data.params,f), continue; end
        MC.(f) = PARAMS{j,3}(data.params.(f));
    end
    
    if isempty(data.data), return; end
    
    % Get timestep from first period tag
    tt = textscan(data.data{1}{1},'%f-%f-%fT%f:%f:%f/%f-%f-%fT%f:%f:%f','CollectOutput',true);
    dt = round(diff(datenum(reshape([tt{:}],6,2)'))*1440);
    UNIT_SCALE = UNIT_SCALE(dt);
    
    if isfield(MC,'dt')
        assert(isequal(MC.dt,minutes(dt)),'Inconsistent timestamps & summarization'); 
    else
        MC.dt = minutes(dt);
    end
    MC.t = (MC.date_begin:MC.dt:MC.date_end-MC.dt/2)' + MC.dt/2;
    MC.interval = 'c';
    
    %MC.t = (datenum(MC.date_begin):dt/1440:datenum(MC.date_end)-dt/2800)' + dt/2880;
    assert(numel(MC.t) == numel(data.data{1}),'Something went wrong with time-labels');
    
    ALL_VARS = matlab.lang.makeValidName(data.headers);
    allvars = cellfun(@lower,ALL_VARS,'unif',0);
    
    reqvars = cellfun(@lower,VAR_LIST_IN,'unif',0);
    [known,col_idx] = ismember(reqvars,allvars);
    
    assert(all(known | OPTIONAL),'Unexpected header! (missing %s) update wps_CAMS_McClear.m',...
        shortliststr(VAR_LIST_IN(~(known|OPTIONAL)),'field',3));

    for j = 1:numel(known)
        if ~known(j) || isempty(VAR_LIST_OUT{j}), continue; end
        MC.(VAR_LIST_OUT{j}) = data.data{col_idx(j)}*UNIT_SCALE(j);
    end
    
    function t = TZdate(x)
        if isfield(MC,'time_ref'), TZ = MC.time_ref; end % (§)
        t = parsetime(x,'InputFormat','yyyy-mm-ddTHH:MM:SS','TimeZone',TZ);
        if isempty(TZ), TZ = t.TimeZone; end
    end
end
    
function dt = parsetimestep(x)
% Parse strings of the form: '0 year 0 month 0 day 0 h 15 min 0 s'
    
    try
        n = regexp(x,'(\d+) year (\d+) month (\d+) day (\d+) h (\d+) min (\d+) s','tokens');
        n = cellfun(@str2double,[n{:}]);
        dt = hours([8765.82, 8765.82/12, 24, 1, 1/60, 1/3600])*n';
    catch
        dt = hours(NaN);
    end
end

function checkfile(file)
    F  = dir(file);
    if isempty(F), throwAsCaller(MException('readcams:nofile','File not found')); end
    if F.bytes == 0, throwAsCaller(MException('readcams:empty','File is empty')); end
    if F.bytes < 10e3
        s = fileread(file);
        code = regexpi(s,'<.*exceptionCode="(.*?)"','tokens');
        if ~isempty(code)
            msg = regexpi(s,'ExceptionText>(.*?)<','tokens');
            msg = [code{1}{1} ': ' msg{1}{1}];
            throwAsCaller(MException('readcams:errorfile',msg));
        end
    end
end

