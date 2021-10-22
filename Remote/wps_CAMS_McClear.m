function [MC,filename,cmd_wps] = wps_CAMS_McClear(Loc,interval,dt,varargin)
% MC = WPS_CAMS_MCCLEAR(LOC,INTERVAL,DT) - Use McClear WPS to retrieve clear-sky
%   irradiation data from http://www.soda-pro.com/web-services/radiation/cams-mcclear
%
%   It is recommended to use WPS_CAMS_MCCLEAR through GETREMOTE('MCCLEAR',..), which makes a local
%   file search first to avoid duplicate requests, and provides more advanced argument options.
%
% INPUT:
%   LOC - Location structure or vector {latitude (DEG),longitude (DEG),altitude (MASL)}
%   INTERVAL - 2 element date-vector (passed through PARSETIME). See note (*)
%   DT - scalar DURATION or time step in minutes (1, 15, 60 or 1440) - also recognized are
%       'h' (hour) and 'd' (day
%
%   (*) NOTE: Data can only be requested as whole-UTC-day chunks. Interval will be converted to
%   UTC unless 'time','TST' is used (see below).
%
% MC = WPS_CAMS_MCCLEAR(...,'email','name@domain.com') - set user of www.soda-pro.com
% MC = WPS_CAMS_MCCLEAR(...,'time','UT') - time reference 'UT' or 'TST'
% MC = WPS_CAMS_MCCLEAR(...,'server',MIROR) - SoDa server. Default is 'http://www.soda-is.com'.
% MC = WPS_CAMS_MCCLEAR(...,'full',true) - retrieve aerosols, albedo, water content, etc.
%   only applicable to 'time' = 'UT' and 1-min resolution. 'full' mode defaults to true for
%   these conditions (verbose mode on CAMS service)
%
% [..] = WPS_MERRA2(..,'-verbose') - provide progress info messages (not the same as '-full')
% [..] = WPS_MERRA2(..,'MaxDist',D) - change default distance D (4.5e4 m) used for validation.
% [..] = WPS_MERRA2(..,'filename',F) - provide file name, otherwise automatically generated from
%   using METEOFILENAME (e.g. mcclear_21N08_104W61_PT60M_2018.csv)
%
% OUTPUT:
%   MC : header structure, result of READCAMS(F,'data',false)
%   F : file name
%   REQ: WPS request
%
% See also: GETREMOTE, READ_CAMS, WPS_MERRA2, COMPLETEMETEODATA

    % Updated for version 3.1, 5.02.2019
    WPS_REQ{1} = '/service/wps?Service=WPS&Request=Execute&Identifier=get_mcclear&version=1.0.0&RawDataOutput=irradiation&DataInputs=';
    WPS_REQ{2} = 'latitude=%.6f;longitude=%.6f;altitude=%.1f;date_begin=%s;date_end=%s;time_ref=%s;summarization=%s;username=%s';

    narginchk(4,Inf);
    parselocation(Loc,'optional',{'TimeZone','name'});

    % parse time-step (*)
    if isa(dt,'duration'), dt = minutes(dt); end
    switch(dt)
        case 1, dt_str='PT01M'; % min
        case 15, dt_str='PT15M';
        case {'h', 60}, dt = 60; dt_str='PT01H';
        case {'d',1440}, dt = 1440; dt_str='P01D';
        otherwise
            error('Incorrect time-step parameter');
    end
    dt = minutes(dt);

    [opt,varargin] = getflagoptions(varargin,{'-verbose','-full'});
    opt.email = '';
    opt.full = opt.full || minutes(dt) == 1;
    opt.time = 'UT';
    opt.server = 'http://www.soda-is.com';
    opt.MaxDist = 1000; % should account for any rounding errors in lat, lon to 0.01°
    opt.filename = meteofilename(Loc,interval(1),interval(2),dt,'prefix','mcclear','suffix','.csv'); 
    opt = getpairedoptions(varargin,opt,'restchk'); 
    
    if isempty(opt.email)
        try
            opt.email = useraccounts('mcclear');
        catch
            error('Unable to retrieve valid credentials'); 
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
    
    str_email = strrep(opt.email, '@', '%2540');

    % McClear WPS expects time reference = UT or TST
    switch upper(opt.time)
        case 'TST',opt.time = 'TST'; TZ = '';
        case {'TU','UT','UTC'}, opt.time = 'UT'; TZ = 'UTC';
        otherwise
            error('Unknown time reference');
    end
    interval = parsetime(interval,'TimeZone',TZ);
    assert(numel(interval) == 2 && diff(interval) >= 0,'Bad date format(s)');

    if opt.full && (dt ~= 1 || ~strcmp(opt.time,'UT'))
        warning('Full mode requires 1-min time-steps at UT reference. Turning it off.');
        opt.full = false;
    end

    verboseprint = @(varargin) opt.verbose && fprintf(varargin{:});
    verboseprint('Attempting to retrieve CAMS McClear data from WPS...\n');

    % Build McClear WPS request
    cmd_wps = [opt.server,WPS_REQ{1},WPS_REQ{2}];
    cmd_wps = sprintf(cmd_wps,Loc.latitude, Loc.longitude, Loc.altitude, ...
                    datestr(interval(1),'yyyy-mm-dd'), datestr(interval(2),'yyyy-mm-dd'),...
                    opt.time, dt_str,str_email);

    if opt.full
        cmd_wps = sprintf('%s;verbose=%d', cmd_wps, opt.full);
    end

    if ispc
        cmd_wget = sprintf('%s -O "%s" "%s"', which('wget.exe'),filename, cmd_wps);
    else %if isunix
        cmd_wget = sprintf('wget -O "%s" "%s"', filename, cmd_wps);
    end
    %cmd_wget = ['wget -O ' filename ' -nv "' cmd_wps '"'];
    [status, wget_res] = system(cmd_wget);
    assert(status == 0,'wget failed with error message: %s\n',wget_res);

    verboseprint('Data saved to file: %s\n',relativepath(filename));

    % Read file header and check that it matches request
    try
        MC = read_CAMS(filename,'data',false);
    catch ERR
        if strcmp(ERR.identifier,'readcans:errorfile'), delete(filename); end
        throw(ERR);
    end

    d = solarposition.arcdist(MC.latitude,MC.longitude,Loc.latitude,Loc.longitude,6.4e6);
    if d > opt.MaxDist
        warning('Location in retrieved file doesn''t match request!');
    end
    if MC.summarization ~= dt
        warning('Time-step in retrieved file doesn''t match request!');
    end
    if MC.date_begin > interval(1) || MC.date_end < interval(2)
        warning('Interval in retrieved file doesn''t match request!');
    end   
    
    verboseprint('Retrieved file seems in order\n');
end
