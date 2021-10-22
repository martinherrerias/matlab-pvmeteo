function [MD,filename,cmd_wps] = wps_MERRA2(Loc,interval,dt,varargin)
% [MD,F,REQ] = WPS_MERRA2(LOC,INTERVAL,DT) - Use CAMS WPS to retrieve MERRA2 weather data.
%   See: http://www.soda-pro.com/web-services/meteo-data/merra
%
%   It is recommended to use WPS_MERRA2 through GETREMOTE('MERRA2',..), which makes a local
%   file search first to avoid duplicate requests, and provides more advanced argument options.
%
% INPUT:
%   LOC - Location structure or vector {latitude (DEG),longitude (DEG),altitude (MASL)}
%   INTERVAL - 2 element date-vector (passed through PARSETIME). See note (*)
%   DT - scalar DURATION or time step in minutes (1, 15, 60 or 1440) - also recognized are
%       'h' (hour), 'd' (day), or 'm' (month)
%
% [..] = WPS_MERRA2(..,'-verbose') - supress all output.
% [..] = WPS_MERRA2(..,'MaxDist',D) - change default distance D (4.5e4 m) used for validation.
% [..] = WPS_MERRA2(..,'filename',F) - provide file name, otherwise automatically generated from
%   using METEOFILENAME (e.g. merra2_21N08_104W61_PT60M_2018.csv)
%
% OUTPUT:
%   MD : header structure, result of READ_MERRA2(F,'data',false)
%   F : file name
%   REQ: WPS request
%
% See also: GETREMOTE, READ_MERRA2, WPS_CAMS_MCCLEAR, COMPLETEMETEODATA

    % WPS_REQ{1} = '/service/wps?Service=WPS&Request=Execute&Identifier=get_mcclear&version=1.0.0&RawDataOutput=irradiation&DataInputs=';
    % WPS_REQ{2} = 'latitude=%.6f;longitude=%.6f;altitude=%.1f;date_begin=%s;date_end=%s;time_ref=%s;summarization=%s;username=%s';

    % Updated for version 3.1, 5.02.2019
    WPS_REQ{1} = 'http://www.soda-is.com/pub/merra2_request.php?';
    WPS_REQ{2} = 'geopoint=%.6f,%.6f&firstday=%s&lastday=%s&duration=%s';

    narginchk(3,Inf);
    parselocation(Loc,'optional',{'TimeZone','name'});
    
    % parse time-step (*)
    if isa(dt,'duration'), dt = minutes(dt); end
    switch(dt)
        case {1,5,10,15,30,'1','5','10','15','30'}, dt_str=num2str(dt); % min
        case {'h','hour', 60}, dt = 60; dt_str='h';
        case {'d','day', 1440}, dt = 1440; dt_str='d';
        case {'month','m',43829.1},dt = 43829.1; dt_str = 'm';
        otherwise
            error('Incorrect time-step parameter');
    end
    dt = minutes(dt);
    
    interval = parsetime(interval);
    assert(numel(interval) == 2 && diff(interval) >= 0,'Bad date format(s)');

    [opt,varargin] = getflagoptions(varargin,{'-verbose'});
    opt.user = '';
    opt.passwd = '';
    opt.MaxDist = 4.5e4; % MERRA-2 grid is 5/8° longitude by 1/2° latitude
    opt.filename = meteofilename(Loc,interval(1),interval(2),dt,'prefix','merra2','suffix','.csv');
    opt = getpairedoptions(varargin,opt,'restchk');
    
    if isempty(opt.user) || isempty(opt.passwd)
        try
           [opt.user,opt.passwd] = useraccounts('merra2');
           parsestruct(opt,{'user','passwd'},'class','char','nonempty');
        catch
           error('Unable to get valid credentials');
        end
    end
    filename = opt.filename;
    
    assert(ischar(filename) || isstring(filename),'Bad file name');
    if isfile(filename), backupdelete(filename,'-warning'); end
    try
        fID = fopen(filename,'w'); fclose(fID);
    catch ERR
        error('Could not create file for writing: %s',getReport(ERR));
    end
    
    verboseprint = @(varargin) opt.verbose && fprintf(varargin{:});
    verboseprint('Attempting to retrieve MERRA2 data from WPS...\n');

    % Build MERRA2 WPS request
    cmd_wps = [WPS_REQ{1},WPS_REQ{2}];
    cmd_wps = sprintf(cmd_wps,Loc.latitude, Loc.longitude,...
        datestr(interval(1),'yyyy-mm-dd'),datestr(interval(2),'yyyy-mm-dd'), dt_str);

    if ispc
        cmd_wget = [which('wget.exe') ' -O ' filename ' --header="soda-user: ' opt.user ...
            '" --header="soda-passwd: ' opt.passwd '" "' cmd_wps '"'];
    else %if isunix
        cmd_wget = ['wget -O ' filename ' --header="soda-user: ' opt.user ...
            '" --header="soda-passwd: ' opt.passwd '" "' cmd_wps '"'];
    end

    [status, wget_res] = system(cmd_wget);
    assert(status == 0,'wget failed with error message: %s\n',wget_res);

    verboseprint('Data saved to file: %s\n',relativepath(filename));
    
    % Read file header and check that it matches request
    try
        MD = read_MERRA2(filename,'data',false);
    catch ERR
        if strcmp(ERR.identifier,'readmerra2:errorfile'), delete(filename); end
        throw(ERR);
    end
    
    d = solarposition.arcdist(MD.latitude,MD.longitude,Loc.latitude,Loc.longitude,6.4e6);
    if d > opt.MaxDist
        warning('Location in retrieved file doesn''t match request!');
    end
    if MD.summarization ~= dt
        warning('Time-step in retrieved file doesn''t match request!');
    end
    if MD.date_begin > interval(1) || MD.date_end < interval(2)
        warning('Interval in retrieved file doesn''t match request!');
    end
    
    verboseprint('Retrieved file seems in order\n');
end

