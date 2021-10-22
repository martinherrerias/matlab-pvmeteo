function MD = read_MERRA2(file,varargin)
% S = READ_MERRA2(FILE) - Read a MERRA-2 meteorological data file:
% 
%     # Modern-Era Retrospective Analysis for Research and Applications (MERRA), version 2
%     # ...
%     # Site latitude (positive means North);48.695
%     # Site longitude (positive means East);9.196
%     # Site altitude (m);545
%     # Date beginning;2019-04-01
%     # Date end;2020-11-30
%     # ...
%     # Time reference (hour);UT
%     # Summarization (period of integration);Hour (h)
%
%     # Columns:
%     # Temperature (K);Temperature at 2 m above ground
%     # Relative humidity (%);Relative humidity at 2 m above ground
%     # Pressure (hPa);Pressure at ground level
%     # Wind speed (m/s);Wind speed at 10 m above ground
%     # ...
% 
%     # MERRA-2 meteorological data
%     # Date;UT time;Temperature;Relative Humidity;Pressure;Wind speed;...
%     2017-01-01;00:10;267.62;89.07;970.58;2.69;...
%     2017-01-01;00:20;267.55;88.95;970.46;2.71;...
%
% OUTPUT: Note change of units in some variables
% 
%     S.info: cellstring, complete header without # markers.
%     S.latitude: degrees north
%     S.longitude: degrees east
%     S.altitude: mASL
%     S.interval: 'c' / 'month'*
%     S.date_begin: datetime scalar (beginning of first day with data)
%     S.date_end: datetime scalar (end of last day with data)
%     S.dt: (duration scalar*) timestep resolution
%     S.t: [N×1] datenum vector of UTC time-steps (center of interval)*
%     S.Ta: [N×1] (°C) Temperature at 2 m above ground
%     S.RH: [N×1] (%) Relative humidity at 2 m above ground
%     S.Patm: [N×1] (Pa) Pressure at ground level
%     S.vw: [N×1] (m/s) Wind speed at 1 m above ground (!) Roughness Length 0.5 m
%     [S.GHI: [N×1] (W/m²) Global Horizontal Irradiance (!)] - not in all versions!
%     S.(..) - {windir,rainfall,snowfall,snowdepth} == {WindDirection, Rainfall, Snowfall, SnowDepth} 
%       [N×1] vectors with unchanged units.
%
%   (*) NOTE: files with monthly averages will result in S.interval = 'month', a vector of
%     S.dt with month durations for each time-step, and time labels S.t at the representative
%     month days SOLARPOSITION.MDAYS(..) for each month.
%
% S = READ_MERRA2(FILE,VARARGIN) - pass additional arguments to READTXTFILE, e.g....
% S = READ_MERRA2(FILE,'data',false) - reads and returns only header-data information.
%
% REF: http://www.soda-pro.com/web-services/meteo-data/merra
%      http://gmao.gsfc.nasa.gov/reanalysis/MERRA-2
% 
% SEE ALSO: WPS_MERRA2, WPS_CAMS_MCCLEAR, COMPLETEMETEODATA

    if nargin < 1, file = '*'; end
    file = pickfile(file);

    WSF = log(1/0.5)/log(10/0.5);  % Wind shear factor, 1 m above ground Roughness Length 0.5 m (!)
    utcdate = @(x) datetime(x,'TimeZone','UTC');

    COLNAMES = {'Date','UTTime','Temperature','RelativeHumidity','Pressure','WindSpeed','WindDirection','Rainfall','Snowfall','SnowDepth','Short_waveIrradiation'};
    FLDNAMES = {''    ,''      ,'Ta'         ,'RH'              ,'Patm'    ,'vw'       ,'windir','rainfall','snowfall','snowdepth','GHI'};
    UNITCHF =  {[]    ,[]      ,@(x) x-273.15,@(x) x/100        ,@(x) x*100,@(x) x*WSF ,[]      ,[]        ,[]        ,[]         , @Wh2W};

    % # PARAMNAME;VALUE pairs in file header
    PARAMS(:,1) = {'.*latitude.*','.*longitude.*','.*altitude.*','.*summarization.*','.*date.?begin.*','.*date.?end.*'};
    PARAMS(:,2) = {'latitude','longitude','altitude','dt','date_begin','date_end'}';
    PARAMS(:,3) = {@str2double,@str2double,@str2double,@parsedt,utcdate,utcdate}';
    
    checkfile(file); % Check if file is not there, empty, or just an error line

    % Read actual file
    data = readtxtfile(file,PARAMS(:,1:2),'regexpmatch',true,'colon',';','ignorecase',true,varargin{:});
    assert(~isempty(data),'Unable to read MERRA file');
    
    data.headers = matlab.lang.makeValidName(data.headers);
    assert(isempty(setdiff(data.headers,COLNAMES)),'Failed to find header: wrong column names?');
    
    MD.info = data.comments;
    
    % Retrieve PARAMNAMES from data.params
    missing = {};
    for j = 1:size(PARAMS,1)
        if ~isfield(data.params,PARAMS(j,2)), missing = cat(1,missing,PARAMS(j,2));
        else
            MD.(PARAMS{j,2}) = PARAMS{j,3}(data.params.(PARAMS{j,2}));
        end
    end
    assert(isempty(missing),'Could''nt find %s',shortliststr(missing,'parameter'));
    
    MD.date_end = MD.date_end + days(1); % MERRA2 uses start of last day!

    if isempty(data.data), return; end % e.g. if (..,'data',false) to get header only
    
    Nt = numel(data.data{1});
    MD.t = parsetime([char(data.data{1}),repmat(' ',Nt,1),char(data.data{2})],'TimeZone','UTC');
    
    if days(MD.dt) > 1 
    % Use representative month days
        MD.t = dateshift(MD.t,'start','year') + days(solarposition.mdoy(month(MD.t))-1);
        dt = dateshift(MD.t,'end','month') - dateshift(MD.t,'start','month') + 1;
        
        MD.interval = 'month';
    else
        [MD.t,dt] = parsetime(MD.t,'-regular','-unique','interval','e');
        % dt = unique(datevec(diff(MD.t)),'rows')*[1440*365,1440*30,1440,60,1,1/60]'; % minutes
        % assert(isscalar(dt) && mod(dt,1) == 0 && dt <= 60);
        MD.t = MD.t - dt/2; % time-steps in MERRA-2 are referred to end of interval

        sameday = @(x,y) diff(dateshift([x,y],'start','day')) == 0;
        if ~sameday(MD.t(1),MD.date_begin) || ~sameday(MD.t(end)+dt,MD.date_end) || ...
                ~isequal(dt,MD.dt)
            warning('File header date_begin, date_end & summarization do not match time-steps');
        end
        
        MD.interval = 'c';
    end
    MD.dt = dt;

    % Assign columns to structure fields
    for k = 1:numel(data.headers)
        j = find(strcmpi(data.headers{k},COLNAMES));
        if isempty(FLDNAMES{j}),continue; end
        MD.(FLDNAMES{j}) = data.data{k};

        missing = MD.(FLDNAMES{j}) == -999;
        MD.(FLDNAMES{j})(missing) = NaN;

        % Change units, if required
        if ~isempty(UNITCHF{j}), MD.(FLDNAMES{j}) = UNITCHF{j}(MD.(FLDNAMES{j})); end
    end

    function y = Wh2W(x)
        y = x./hours(dt);
    end
end

function dt = parsedt(x)
    switch x
    case 'Month', dt = years(1)/12;
    case 'Day (d)', dt = days(1);
    case 'Hour (h)', dt = hours(1);
    otherwise
        try
            n = regexp(x,'(\d+) minutes','tokens');
            n = str2double(n{1}{1});
            dt = minutes(n);
        catch
            dt = minutes(nan);
        end
    end
end

function checkfile(file)
    F  = dir(file);
    if isempty(F), throwAsCaller(MException('readmerra2:nofile','File not found')); end
    if F.bytes == 0, throwAsCaller(MException('readmerra2:empty','File is empty')); end
    if F.bytes < 10e3
        s = fileread(file);
        s = strsplit(strtrim(s),newline());
        if isscalar(s) && contains(s{1},'error','ignorecase',true)
           throwAsCaller(MException('readmerra2:errorfile','File contains: %s',s));
        end
    end
end


