function [MD, Time, Location] = getMeteoData(filename,filetype)
% [MD, Time, Location] = GETMETEODATA(filename,[filetype])
% Import Meteorological Data and site Location from file 'filename', of type 'filetype'
%
% filename	(including path) points to the file with the meteo data (Use [] to search manually)
% filetype  ('MeteoNorm_solar','MeteoNorm_min','SolarGIS') specifies input file format. 
%			Use [] to detect format automatically.
%
% 'MeteoNorm_min':
%     Data-Set Title
%     [LatN] ~ [LonE] ~ [Alt] ~ [TZ] ~ [Offset] ~
%     m dm dy h mi G_Gh hs G_Gex G_Gh_hr G_Dh G_Gk G_Dk G_Bn Ta FF G_Gcs
%     [~ ~ DoY H M  GHI  ~   ~   ~     DHI   ~    ~   BNI Ta vw]
%
% 'MeteoNorm_solar':
%     "Data-Set Title"
%     [LatN],[LonE],[Alt],[TZ],[Offset],1,1,"SYN"
%     ~  
%     ~(Monthly-averages)
%     ~(Monthly-averages)
%     ~(Monthly-averages)
%     [GHI DHI Ta windspeed]
%
% 'SolarGIS':
%    #..
%    #Latitude: [LatN]
%    #Longitude: [LonE]
%    #Altitude: [Alt]
%    #UTCoffset: [TZ]
%    #interval: ['mid'='instant'] 'start'/'end'
%    #..
%    GHI DHI BNI Ta...
%    ['dd.mm.yy' 'HH:MM' GHI BNI DHI Ta vw]
%
% 'CAMS': (see READ_CAMS)
%
% 'generic' (as written by WRITEMETEOFILE)
%
%    #..
%    # Name: Example [LOCATION.name]
%    # Latitude: 21.000000 [LOCATION.latitude]
%    # Longitude: -110.000000 [LOCATION.longitude]
%    # Altitude: 1800.0 [LOCATION.altitude]
%    # TimeZone: [LOC.TimeZone, or TIMESTAMPS.TimeZone, if available] (*)
%    #
%    # period: 2020-01-01T00:00:00Z/2020-12-31T23:55:00Z  [ISO 8601] 
%    # interval: yyyy-mm-ddTHH:MM:SSZ/PT5M           [ISO 8601 format] / {point, start, end, center}
%    #
%    # [COMMENTS]
%    #
%    # time GHI DHI BNI Ta...
%    2019-03-11T15:30:00Z 100.0 200.0 50.0 20.0 ..
%    2019-03-11T15:35:00Z 80.0 190.0 30.0 20.0 ..
%
% OUTPUT (additional/missing fields based on file type):
%     Location.name - ID string
%     Location.latitude - scalar, + North, degrees
%     Location.longitude - scalar, + East, degrees
%     Location.altitude - scalar, meters above sea level
% 
%     MD.info - all text in original files
%     MD.timestep - DURATION of summarization interval
%     MD.DoY - Day of Year, from 1.0 to 366.0
%     MD.GHI - Global Horizontal Irradiance
%     MD.BNI - Beam Normal Irradiance
%     MD.DHI - Diffuse Horizontal Irradiance
%     MD.Ta - Ambient temperatue
%     MD.vw - wind speed (if not available returns scalar 1m/s)
% 
%     Time - DATETIME vector.
% 
% See also: COMPLETEMETEODATA, WRITEMETEOFILE

	% UI file search if filename is not provideddata
    if nargin < 1 || isempty(filename), filename = '*'; end
    filename = pickfile(filename,'Select a Meteo-Data file');

    % Attempt to read a file with #comments followed by delimited data
    % Try finding the following #Parameter:value pairs in the comments
    knownparams = {'latitude','latitude';
                   'longitude','longitude';
                   'altitude','altitude';
                   'name','name';
                   'title','name';
                   'time.?(reference|zone)','timezone';
                   'UTCoffset','utcoffset';
                   '(summarization|interval)$','interval'};
    data = readtxtfile(filename,knownparams,'ignorecase',true,'regexpmatch',true,'TreatAsEmpty','NaN');

    if isempty(data) 
    % If it didn't work, import data blindly (just for type recognition)
        data = importdata(filename);
        
        % clear empty rows (default for importdata until R2016b) 
        data.textdata = data.textdata(~all(cellfun(@isempty,data.textdata),2),:);
        
    elseif ~isempty(data.headers)
        data.headers = matlab.lang.makeValidName(data.headers);
    end
    
	% If filetype was not explicitly specified, try to guess
    if nargin<2 || isempty(filetype)
		filetype = guessdatatype(data);
		if filetype == 0
			error('getMeteoData: File type not recognized, or not yet supported');
		end
    end
	
	% Go to importing subroutines based on file type
    switch lower(filetype)
	   %case 'TMY3',                       [MD, Time, Location] = ReadTMY3(filename);
        case {'generic','writemeteofile'}, [MD,Time,Location] = readgeneric(data);
		case {'solargis','gislike'},       [MD, Time, Location] = ReadSolarGIS(data);
        case 'cams'
            data = read_CAMS(filename);
            %Time = parsetime(data.t,'timezone',data.time_ref);
            Time = data.t;
            Location.latitude = data.latitude;
            Location.longitude = data.longitude;
            Location.altitude = data.altitude;
            MD = rmfield(data,{'latitude','longitude','altitude','t','time_ref'});
		case 'meteonorm_solar',         [MD, Time, Location] = ReadMeteoNormSolar(data); 
		case 'meteonorm_min',           [MD, Time, Location] = ReadMeteoNormMin(data);
        otherwise
            error('Unknown file type');
    end
    parselocation(Location);
    
    src = ['Source: ' absolutepath(filename)];
    MD.info = [{src};MD.info];
    
    for f = fieldnames(MD)'
        if isnumeric(MD.(f{1})), MD.(f{1}) = single(MD.(f{1})); end
    end
end

function type = guessdatatype(data)
% tries to find out which format of Meteo-Data file we're dealing with, based on things like
% header lines, number of columns, comments, etc...

    type{1} = 0;
    foundit = false;
    %if ~foundit && checkTMY3(data), type{1} = 'TMY3'; foundit = true; end
    if ~foundit && checkMeteoNormSolar(), type{1} = 'MeteoNorm_solar'; foundit = true; end
    if ~foundit && checkMeteoNormMin(), type{1} = 'MeteoNorm_min'; foundit = true; end
    if ~foundit && checkSolarGIS(), type{1} = 'SolarGIS'; foundit = true; end
    if ~foundit && checkCAMS(),type{1} = 'CAMS'; foundit = true; end
    if ~foundit && checkgeneric(), type{1} = 'Generic'; foundit = true; end

    type = type{1};
    if foundit, fprintf('Attempting to read meteo-file as %s...\n',type); end

    return;

    % function ok = checkTMY3(data)
    % % importdata(data) gets stuck before the header line, so...
    % 	ok = size(data.data,1)==1 && size(data.textdata == 1);  
    % end
    function ok = checkgeneric()
        ok = false;
        % Result must be full-output of readtxtfile()
        if ~isstruct(data) || ~all(isfield(data,{'comments','headers','data'})) || ...
            isempty(data.headers) || isempty(data.data) || isempty(data.comments), return; end
        
        if numel(data.headers) > 1 && contains(data.headers{1},'time','IgnoreCase',true) && ...
            all(isfield(data.params,{'name','latitude','longitude','altitude','interval'}))
            ok = true;
        end
    end
    function ok = checkSolarGIS()
        ok = false;
        % Result must be full-output of readtxtfile()
        if ~isstruct(data) || ~all(isfield(data,{'comments','headers','data'})) || ...
            isempty(data.headers) || isempty(data.data) || isempty(data.comments), return; end
        
        if numel(data.headers) > 5 && ...
            strcmpi(data.headers{1},'date')&&...
            strcmpi(data.headers{2},'time')&&...
            strcmpi(data.headers{3},'GHI')&&...
            any(strcmpi(data.headers{4},{'DNI','BNI'}))&&...
            any(strcmpi(data.headers{5},{'DIF','DHI'}))&&...
            any(strcmpi(data.headers{6},{'TEMP','Ta'}))
            % if all headings match, go for it
            ok = true;
        end
    end
    function ok = checkMeteoNormSolar()
        ok = false;
        % Result must be output of importdata()
        if ~isstruct(data) || ~all(isfield(data,{'textdata','data'})), return; end
        % Data in columns: [GHI DHI Ta vw]
        if ~isstruct(data); return; end
        if data.textdata{1}(1) == '"' && numel(data.textdata) == 5 &&...
            size(data.data,2) >= 4 &&...				% {GHI DHI Ta vw}
            all(data.data(:,1)>=0) && ...				% GHI >= 0
            all(data.data(:,2)>=0) && ...				% DHI >= 0
            all(data.data(:,2)<=data.data(:,1)) &&...	% DHI <= GHI
            all(data.data(:,3)<60)&&...					% Ta < 60ºC
            all(data.data(:,4)>=0)						% vw >= 0
            ok = true;
        end
    end
    function ok = checkMeteoNormMin()
        ok = false;
        % Result must be output of importdata()
        if ~isstruct(data) || ~isfield(data,'colheaders'), return; end
        % The following header names exist:
        if ...
            any(strcmp('dy',data.colheaders))&&...
            any(strcmp('h',data.colheaders))&&...
            any(strcmp('mi',data.colheaders))&&...
            any(strcmp('G_Gh',data.colheaders))&&...
            any(strcmp('G_Gex',data.colheaders))&&...
            any(strcmp('G_Dh',data.colheaders))&&...
            any(strcmp('G_Bn',data.colheaders))&&...
            any(strcmp('Ta',data.colheaders))&&...
            any(strcmp('FF',data.colheaders))
            ok = true;
        end
    end
    function ok = checkCAMS()
        ok = false;
        % Result must be full-output of readtxtfile()
        if ~isstruct(data) || ~all(isfield(data,{'comments','headers','data'})) || ...
            isempty(data.headers) || isempty(data.data) || isempty(data.comments), return; end
        
        if isequal(regexpi(data.headers{1},'observation([_\w]?)period'),1) && ...
            any(strcmp('ClearSkyGHI',data.headers)) && ...
            any(strcmp('ClearSkyDHI',data.headers)) && ...
            any(strcmp('ClearSkyBNI',data.headers))
            
            % if all headings match, go for it
            ok = true;
        end
    end
end
	
% function [MD, Time, Location] = ReadTMY3(filename)
% 
% 	fid = fopen(file);
% 	firstline = textscan(fid,'%d %s %d %d %f %f %f',1,'delimiter',',');
% 	headers = textscan(fid,'%s',68,'delimiter',',');
% 	%data = textscan(fid,'%s %s 
% 	fclose(fid);
% end

function [MD, t, Location] = readgeneric(data)
% Read the result of readtxtfile() for a file written by WRITEMETEOFILE
%
%    #..
%    # Name: Example [LOCATION.name]
%    # Latitude: 21.000000 [LOCATION.latitude]
%    # Longitude: -110.000000 [LOCATION.longitude]
%    # Altitude: 1800.0 [LOCATION.altitude]
%    # [TimeZone: 'IANA Time-Zone']
%    #
%    # period: yyyy-mm-ddTHH:MM:SSZ/yyyy-mm-ddTHH:MM:SSZ  [ISO 8601] 
%    # interval: yyyy-mm-ddTHH:MM:SSZ/PTXX           [ISO 8601] 
%    #
%    # [COMMENTS]
%    #
%    # time GHI DHI BNI Ta...
%    2019-03-11T15:30:00Z 100.0 200.0 50.0 20.0 ..

    % Parse variable names to case-sensitive convention
    data.headers = standardnames(data.headers);

    % Read Location details and UTC offset from comment lines
    REQFIELDS = {'name','latitude','longitude','altitude','interval'};
    missing = setdiff(REQFIELDS,fieldnames(data.params));
    assert(isempty(missing),'Missing %s',shortliststr(missing,'parameter','colon',':'));
    
    REQVARS = {'time'}; 
    missing = setdiff(REQVARS,data.headers);
    assert(isempty(missing),'Missing %s',shortliststr(missing,'column','colon',':'));

    if isfield(data.params,'name'), Location.name = data.params.name; end
    Location.latitude = str2double(data.params.latitude);
    Location.longitude = str2double(data.params.longitude);
    Location.altitude = str2double(data.params.altitude); 
    
    % Remove (if any) brackets that specify non-printed sections of interval string
    data.params.interval = strrep(data.params.interval,'[','');
    data.params.interval = strrep(data.params.interval,']','');

    % Test for ISO 8601 end-summarization format string: e.g. PT5M/yyyy-mm-ddThh:MM:SSZ
    fmt = regexpi(data.params.interval,'^(?<period>PT[HMS\d]*)/(?<fmt>[ymdHMS:T-]+)(?<utc>.*)','names');  % end
    if ~isempty(fmt)
        UTC = fmt.utc;
        period = fmt.period;
        fmt = fmt.fmt;
        data.params.interval = 'e';
    else
    % Test for ISO 8601 beginning-summarization format string: e.g. yyyy-mm-ddThh:MM:SSZ/PT5M
        fmt = regexpi(data.params.interval,'^(?<fmt>[ymdHMS:T-]+)(?<utc>.*)/(?<period>PT[HMS\d]*)','names');
        if ~isempty(fmt)
            UTC = fmt.utc;
            period = fmt.period;
            fmt = fmt.fmt;
            data.params.interval = 'b';
        else
       % Test for ISO 8601 duration+context: e.g. PT5M - centered
            [period,rest] = regexpi(data.params.interval,'(PT[HMS\d]*)','match','split');
            [UTC,rest] = regexpi([rest{:}],'([+-][\d:]{2,5})|Z','match','split');
            [fmt,rest] = regexpi([rest{:}],'([ymdHMS:T-]){5,}','match','split'); 
            rest = regexpi([rest{:}],'[a-z]*','match'); % remove special characters
            data.params.interval = strjoin(rest,' ');
            period = [period{:}];
            UTC = [UTC{:}];
            fmt = [fmt{:}];
        end
    end
    % Parse UTC offset (Z) or (UTC+00:00)
    if ~isempty(UTC), UTC = checktimezone(UTC,'UTC'); end
    
    MD.info = data.comments;

    data.params.interval = checksummarization(data.params.interval);

    % Convert date strings to DATETIME, check regularity, guess dt, and set interval-summarization
    % to center.
    [t,dt] = parsetime(cat(1,data.data{1}{:}),'TimeZone',UTC,'InputFormat',fmt,...
        'interval',{data.params.interval,'c'});
    
    MD.interval = data.params.interval;
    if MD.interval ~= 'i', MD.interval = 'c'; end 
    
    % Parse ISO 8601 interval duration PT[nH][nM][nS] -> minutes
    if ~isempty(period)
        period = regexpi(period,'PT(?<hour>\d+(?=H))?H?(?<min>\d+M)?(?<sec>\d+S)?','names');
        period = str2double({period.hour(1:end-1),period.min(1:end-1),period.sec(1:end-1)});
        period(~isfinite(period)) = 0;
        period = hours(period(1)) + minutes(period(2)) + seconds(period(3));
    else
        period = minutes(NaN); 
    end

    if isfinite(period) && isfinite(dt)
        if dt ~= period
            warning(['Inconsistent time-step duration: %s (header) and %s (diff), ',...
                     'Using %s for interval-labeling offset correction'],...
                     char(period),char(dt),char(min(period,dt)));
        end
    end
    period = min(dt,period); % keep one if the other is NaN
    MD.timestep = period;
    
    if isfield(data.params,'timezone')
        data.params.timezone = checktimezone(data.params.timezone);
        if ~isempty(data.params.timezone)
            Location.TimeZone = data.params.timezone;
        end
    end  

    % ... and just dump the rest under a field named as each header
    for j = 2:numel(data.headers), MD.(data.headers{j}) = data.data{j}; end
end

function [MD, t, Location] = ReadSolarGIS(data)
% Read the result of readtxtfile() for a file with format:
%
%    #..
%    #Latitude: [LatN]
%    #Longitude: [LonE]
%    #Altitude: [Alt]
%    #UTCoffset: [TZ]
%    #..
%    GHI DHI BNI Ta...
%    ['dd.mm.yy' 'HH:MM' GHI BNI DHI Ta vw]

    if ~isfield(data.params,'altitude')
    % Look for altitude in comment: # ...&loc=+00.000000,00.000000&z=00 
        z = regexpi(data.comments,'(?:loc=[+\-\.,0-9]*&z=)(\d*)','tokens');
        while iscell(z) && ~isempty(z), z = cat(1,z{:}); end
        if isempty(z), z = 'NaN'; end
        data.params.altitude = z;
    end
    data.params = completestruct(data.params,...
        struct('latitude','NaN','longitude','NaN','utcoffset','','name','unknown'));
    
    Location.name = data.params.name;
    Location.latitude = str2double(data.params.latitude);
    Location.longitude = str2double(data.params.longitude);
    Location.altitude = str2double(data.params.altitude);

    UTCoffset = checktimezone(data.params.utcoffset,'UTC');
    
    % Parse variable names to case-sensitive convention
    data.headers = standardnames(data.headers);
    
    REQVARS = {'date','time','GHI','BNI','DHI','Ta'}; 
    missing = setdiff(REQVARS,data.headers);
    assert(isempty(missing),'Missing %s',shortliststr(missing,'column','colon',':'));

    [MD.interval,offset] = checksummarization(data.params.interval);
    if MD.interval == 'i', offset = 0; else, offset = offset - 0.5; end

    % Find the number of the header line
    MD.info = data.comments;

    % convert date strings to numbers (Day count)
    Nt = size(data.data{1},1);
    [t,dt] = parsetime([cat(1,data.data{1}{:}),repmat(' ',Nt,1),cat(1,data.data{2}{:})],...
        'InputFormat','dd.mm.yyyy HH:MM','TimeZone',UTCoffset);
    if offset~=0
        t = t - offset.*dt; % set interval summarization to center
    end
    MD.timestep = dt;

    % ... and just dump the rest under a field named as each header
    for j = 3:numel(data.headers), MD.(data.headers{j}) = data.data{j}; end

    % if ~isfield(MD,'vw')
    %     MD.vw = [];
    %     warning('getMeteoData:noWindSpeed','Windspeed data not provided')
    % end
end

function [MD, t, Location] = ReadMeteoNormMin(data)
% Reads a MeteoNorm 7.0 'Standard minute' file with format:
%
%   Data-Set Title
%   [LatN] ~ [LonE] ~ [Alt] ~ [TZ] ~ [Offset] ~
%   m dm dy h mi G_Gh hs G_Gex G_Gh_hr G_Dh G_Gk G_Dk G_Bn Ta FF G_Gcs
%   [~ ~ DoY H M  GHI  ~  ~   ~     DHI   ~    ~   BNI Ta vw]
% 
% Where [x] stands for numeric fields. [LatN], [LonE], and [Alt] (latitude/longitude/altitude) will 
% be stored directly in Location (must match pvlmod_spa conventions), [TZ] (time zone) is passed along 
% with Time structure. For HOURLY-DATA-ONLY (!) [Offset] (minutes) are applied to all time-stamps  
% in the file.
% (!) One-minute files also show a -30 min offset which seems to be wrong.
% 
% For requirements of Location and UTCoffset read documentation for pvlmod_spa
% and pvl_alt2pres.

    % Read Location from header lines
    Location.name = data.textdata{1};
    LocationData = regexp(data.textdata{2},repmat('([\d.]+)[^\d]+',1,5),'tokens');
    LocationData = str2double(LocationData{1});
        Location.latitude = LocationData(1);
        Location.longitude = LocationData(2);
        Location.altitude = LocationData(3);
        UTCoffset = LocationData(4);
        MinOffset = LocationData(5);
    clear LocationData
    txt = data.textdata; data = data.data;
    
    if ~all(data(:,5) == 0), MinOffset = 0; end  % PROVISIONAL: possible error in MeteoNorm Output

    % Read numeric data according to columns:
    % [~ ~ DoY(int) H M  GHI  ~  ~   ~     DHI   ~    ~   BNI Ta vw ]
    DoY = data(:,3) + data(:,4)/24 + (data(:,5) + MinOffset)/1440;
    t = doy2time(DoY,UTCoffset);

    MD.GHI = data(:,6);
    MD.DHI = data(:,10);
    MD.BNI = data(:,13);
    MD.Ta = data(:,14);
    MD.vw = data(:,15);
    MD.info = txt;
    MD.timestep = t(2)-t(1);
end

function [MD, t, Location] = ReadMeteoNormSolar(alldata)
% Reads a MeteoNorm 7.0 'PVsyst' file with the follwing format:
%
%     "Data-Set Title"
%     [LatN],[LonE],[Alt],[TZ],[Offset],1,1,"SYN"
%     ~  
%     ~(Monthly-averages)
%     ~(Monthly-averages)
%     ~(Monthly-averages)
%     [GHI DHI Ta windspeed]
%
% Where [x] stands for numeric fields. [LatN], [LonE], and [Alt] (latitude/longitude/altitude) will 
% be stored directly in Location (must match pvlmod_spa conventions), [TZ] (time zone) is passed along 
% with Time structure. [Offset] is currently ignored.
%
% For requirements of Location and UTCoffset read documentation for pvlmod_spa
% and pvl_alt2pres.

 txt = alldata.textdata; data = alldata.data;

% Read Location from header lines
    Location.name = txt{1}(2:end-1);

    LocationData = sscanf(txt{2},'%f,%f,%f,%d');
        Location.latitude = LocationData(1);
        Location.longitude = LocationData(2);
        Location.altitude = LocationData(3);
        UTCoffset = LocationData(4); % UTC offset sign seems reversed
    clear LocationData
    
    n = size(data,1);
    assert(mod(n,365) == 0,'Expecting time-steps to be integer multiple of 365');
    t = doy2time(1+365*(0:n-1)/n,UTCoffset)';
%     t = doy2time(1+365*(0.5:n)/n,UTCoffset)';
    
    % NOTE: it's not clear what the summarization convention inside PVsyst is. Using centered
    %   intervals causes UTC errors of {0,1,-1} hour + 30 min in most cases. Start of interval
    %   (above) results in detected integer UTC errors of {0,1,-1} hour.
    
    MD.GHI = data(:,1);
    MD.DHI = data(:,2);
    MD.Ta = data(:,3);
    MD.vw = data(:,4);
    MD.info = txt;
    MD.timestep = t(2)-t(1);
end

function t = doy2time(DoY,UTCoffset)
% Convert fractional Day-of-Year values (1.0 = Jan 1. 00:00, 366.0 = Dec 31. 24:00) into
% a DATETIME vector, using the current year (or last, avoiding leap-years).

    UTCoffset = checktimezone(UTCoffset,'UTC'); 

    % Use this year or last (non-leap year)
    y = year(now()-60); % -2 months ensures MERRA2 data is available
    y = y - (rem(y,4) == 0 && (rem(y,100) ~= 0 || rem(y,400) == 0));
    
    t = datetime(datenum(y,1,1) + DoY - 1,'convertfrom','datenum','TimeZone',UTCoffset);
end

function output = standardnames(input)
% Standardize (caps-sensitive) variable names, replacing known aliases, and removing spaces,
% punctuation, and upper-cases from any other names.
%
% FUTURE: consider checking units, transform equivalent fields (e.g. sza & el), use regexp?

    assert(iscellstr(input) || isstring(input));

    % Known aliases for meteorological variables (case-insensitive, left-most is official)
    ALIAS = {
        {'time','timestamp','datetime','UTTime'};
        {'AMa','AM','airmass','airmassabsolute'};
        {'Az','SunAz','azimuth','solarazimuth'};
        {'sza','sza','zenith'};
        {'El','SunEl','elevation'};
        {'CSBHI','ClearSkyBHI'};
        {'CSBNI','ClearSkyBNI','ClearSkyDNI'};
        {'CSDHI','ClearSkyDHI'};
        {'CSGHI','ClearSkyGHI'};
        {'GHI','ShortwaveIrradiation'};
        {'BNI','DNI'};
        {'DHI','DIF'};
        {'ENI','Gextra','DNITOA'};
        {'TOA','TOA','GHITOA'};
        {'Patm','Pressure'};
        {'QcFlag','flag'};
        {'RH','RelativeHumidity','RelHum'};
        {'Ta','Temperature','Tamb','temp','temperatureair'};
        {'Tm','ModuleTemp','ModuleTemperature','temperaturemodule'};
        {'tpw','tcwv','totalprecipitablewater'};
        {'vw','WindSpeed'};
        {'windir','WindDirection'};
        {'kd','fractiondiffuse'};
        {'kt','clearnessghi'};
    };

    % Remove spaces & punctuation
    output = lower(matlab.lang.makeValidName(input));
    output = regexprep(output,'[\s-_]','');

    used = false(numel(ALIAS),1);
    for k = 1:numel(output)
        for j = find(~used)'
            if isempty(j), continue; end
            if any(strcmpi(output{k},ALIAS{j}))
               output(k) = ALIAS{j}(1);
               used(j) = true;
               break;
            end
        end
        if all(used), break; end
        if used(j), continue; end
    end

end
