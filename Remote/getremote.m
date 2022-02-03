function D = getremote(type,Loc,interval,dt,varargin)
% MC = GETREMOTE(TYPE,LOC,INTERVAL,DT) - retrieve McClear/MERRA2 data for a given location LOC 
%   and INTERVAL, and with a given resolution DT.
%
% MC = GETREMOTE(TYPE,LOC,T,DT,'-sync') - retrieve McClear/MERRA2 data for timesteps T. In this
%   case, DT doesn't have to match the time-step resolution of T. It sets the minimum desired 
%   resolution of the TYPE data, that will be interpolated/downsampled as required.
%
%   Beyond being a wrapper for WPS_MERRA2 and WPS_CAMS_MCCLEAR, GETREMOTE is thought as a tool to
%   use and mantain a library of files already downloaded, using them whenever possible, and only
%   launching web requests as a last resource.
%
%   1. Search for existing PATH/TYPE*.csv files in PWD and the current MATLAB PATH (see 'path'
%      and 'pattern' options below).
%   2. Filter out files whose name, based on the standard of METEOFILENAME, does not match the
%      request. That is files for locations further than 'MaxDist', or with time-steps larger 
%      (or not integer divisors) of DT.
%   3. Read the header of short-listed files, and filter again. On UI mode, offer to rename these
%      files according to METEOFILENAME, for fast filterning on the next run.
%   4. On UI mode, provide UI pick of additional files (if required)
%   5. If requirements are still not covered, launch web request, and verify download.
%   6. Merge all of the retrieved files to cover the requested interval, downsampling if necessary,
%      and appling pixel-altitude correction for MERRA2 data.
%   7. Filter / interpolate to timesteps T, if required.
%
% INPUT:
%   TYPE - 'mcclear' or 'merra2', see READ_CAMS and READ_MERRA2 for details
%   LOC - Location structure (see PARSELOCATION).
%   INTERVAL - 2-vector datetime,  date-number or date-vector, parsed by PARSETIME
%   DT - time step DURATION, or number (minutes: 1, 15, 60 or 1440) - also recognized are 
%       'h' (hour), 'd' (day), and for MERRA2 data, 'm' (month)
%
% OUTPUT:
%   MC : MeteoData object. Variable names are the same as in READ_CAMS / READ_MERRA2 and should
%       match the conventions in METEODATA (see METEODATA.WHATIS). Variable sources, however, 
%       are prefixed by 'mcclear.*' / 'merra2.*' so that the resulting object can be merged to
%       other MeteoData objects without conflicts.
%
% Example: 15-min data for last month
%   MC = WPS_CAMS_MCCLEAR([45.0 25.0],now-32,now-2,15)
%
% See also: READ_CAMS, READ_MERRA2, WPS_MERRA2, WPS_CAMS_MCCLEAR, COMPLETEMETEODATA

    narginchk(4,Inf);
    
    TYPES = {'merra2','mcclear'};
    assert(ischar(type) && any(strcmpi(type,TYPES)),'Expecting MERRA2 or MCCLEAR');
    type = TYPES{strcmpi(type,TYPES)};
    
    [opt,varargin] = getflagoptions(varargin,{'-sync','-verbose','-nofix','-ignorenames'});
    opt.UI = runningfromUI();
    opt.path = {};
    opt.pattern = {};
    opt.MaxRequests = 24;
    opt.interval = 'c';
    opt.MaxDist = 1000*(1+44*isequal(type,'merra2'));
    opt.tol = 0; % allow some fraction of the data to be missing
    
    [opt,varargin] = getpairedoptions(varargin,opt);

    % Parse location, try to autocomplete from GEONAMES, skip missing fields
    Loc = parselocation(Loc,'-useweb');
    if ~isfield(Loc,'altitude'), Loc.altitude = -999; end
    if ~isfield(Loc,'TimeZone'), Loc = parselocation(Loc,'-soft'); end
    
    interval = parsetime(interval,'TimeZone','keep');
    if isempty(interval.TimeZone), interval.TimeZone = Loc.TimeZone; end

    if opt.sync 
        t = interval;
        if opt.interval == 'i'
            [t,~,idx,ia] = parsetime(t,'interval','i','-sorted','-unique','gridded',false);
            data_dt = 0;
        else
        % Bring time-labels to center of interval
            [t,data_dt,idx,ia] = parsetime(t,'-regular','interval',{opt.interval,'c'});
            if isempty(dt), dt = data_dt; end
        end
        t.TimeZone = 'UTC';
        interval = [t(1)-dt/2,t(end)+dt/2];
    else
        assert(~isempty(dt) && isnat(dt) && dt > 0, 'Incorrect time-step parameter');
        assert(numel(interval) == 2 && diff(interval) > 0,'Bad interval');
    end
    dt = checkstep(dt,type);

    if ~isduration(opt.tol), opt.tol = diff(interval)*opt.tol; end
    UTCinterval = useful_UTC_interval(interval,Loc); % span of whole UTC days
    
    verboseprint = @(varargin) opt.verbose && fprintf(varargin{:});
    verboseprint('Request for %s: %s\n',upper(type),meteofilename(Loc,UTCinterval(1),UTCinterval(2),dt));
    
    % Attempt to find local files for the specified request...
    opt.justguessing = true;
    files = searchexisting(type,opt.path,opt.pattern);
    [files,cover] = filtercandidates(files,type,Loc,UTCinterval,dt,opt);
    verboseprint([nthings(numel(files),'file') ' in search directories\n']);
    
    opt.ignorenames = true;
    opt.justguessing = false;
    
    missing = cover == 0;
    if opt.UI && any(missing)
    % Complete with UI selected files, if available...
    
        [~,msg] = listgaps(missing,UTCinterval);
        msg = ['Pick ' type ' file(s) (cancel to retrieve from WPS): ',msg];
        newfiles = pickfile('*.csv',Inf,msg,'ui',1,'-fullpath');
        [newfiles,newcover] = filtercandidates(newfiles,type,Loc,UTCinterval,dt,opt);
        if ~isempty(newfiles)
            [files,cover] = mergedcover(files,cover,newfiles,newcover);
            verboseprint(nthings(numel(files),' UI selected file\n'));
            missing = cover == 0;
        end
    end
    
    if any(missing)
    % Finally, retrieve missing data from WPS
        gaps = listgaps(cover == 0,UTCinterval);
        newfiles = remote_requests(type,Loc,gaps,dt,opt,varargin);
        [newfiles,newcover] = filtercandidates(newfiles,type,Loc,UTCinterval,dt,opt);
        if ~isempty(newfiles)
            [files,cover] = mergedcover(files,cover,newfiles,newcover);
            verboseprint([nthings(numel(files),'file') ' from WPS\n']);
            missing = cover == 0;
        end
    end
    % assert(all(cover),'');
    if any(missing)
        [~,msg] = listgaps(missing,UTCinterval);
        assert(nnz(missing)*dt <= opt.tol,msg);
        warning(msg);
        if all(missing)
            D = MeteoData();
            return; 
        end
    end
    
    % Read, [downsample,] merge, and clip data to required interval.
    [D,M] = mergedata(type,files,cover,UTCinterval,dt);
    verboseprint('Using %s:\n',nthings(numel(files),'file','-noN'));
    verboseprint('%s\n',shortliststr(files,'','-newlines'));
    
    % Apply altitude correction
    if strcmpi(type,'merra2') && ~opt.nofix
        if isfield(M,'altitude'), data_h = M.altitude;
        else, data_h = [M.files.altitude]; data_h = data_h(D.file); 
        end
        [D,msg] = altitude_correction(D,data_h,Loc.altitude); 
        if ~isfield(M,'info'), M.info = {}; end
        M.info = cat(1,M.info,msg);
        verboseprint('%s\n',strjoin(msg,newline()));
    end
    
    if opt.sync
       D = sync(D,t,data_dt,idx(ia));
    end
    
    if isempty(D)
        warning('Resulting table contains no useful data');
    end

    D.Properties.UserData = M;
    D = MeteoData(D,'interval','c','location',Loc,'useweb',false);
    D.source = cellfun(@(x) {strjoin([{type},x],'.')},fieldnames(D),'unif',0);
end

function dt = checkstep(dt,type)
% Parse DT (return a scalar DURATION). Make sure it's a valid time-step according to TYPE, or pick
% the next lower valid time-step, with a warning.

    switch type
    case 'merra2'
        VALID = minutes([1,5,10,15,30,60,1440,43829.1]);
        default = minutes(60);
    case 'mcclear'
        VALID = minutes([1,15,60,1440]);
        default = minutes(15);
    end
    
    if isempty(dt)
        dt = default;
        warning('Using default resolution (%s) for %s',isoduration(dt),type);
        return;
    end
    if ~isduration(dt)  
        if isnumeric(dt), dt = minutes(dt);
        elseif ischar(dt)
            switch lower(dt)
            case {'1','5','10','15','30','60'}, dt = minutes(num2str(dt));
            case {'h','hour'}, dt = hours(1);
            case {'d','day'}, dt = days(1);
            case {'month','m'}, dt = years(1)/12;
            otherwise
                try dt = isoduration(dt,true); 
                catch
                    dt = minutes(NaN);
                end
            end
        else, dt = minutes(NaN);
        end
    end
    assert(dt > 0,'Incorrect time-step parameter');
    
    if ~any(dt == VALID)
        if dt < VALID(1), j = 1; else, j = find(dt < VALID,'last'); end
        warning('%s is not a standard time step, using %s',isoduration(dt),isoduration(VALID(j)))
        dt = VALID(j);
    end
end

function files = searchexisting(type,paths,pattern)
% Search for files of the  given TYPE in PWD & project directories.
    
    if isempty(paths)
        paths = strsplit(path(),':');
        paths(~contains(paths,type,'ignorecase',true)) = [];
    end

    if isempty(pattern)
        switch type
        case 'merra2'
            pattern = {'merra2*.csv','MERRA2*.csv','SoDa_MERRA2*.csv'};
        case 'mcclear'
            pattern = {'McClear*.csv','mcclear*.csv','MCCLEAR*.csv','irradiation-*.csv'};
        end
    end
    
    paths{end+1} = pwd();
    if any(strcmp(who('global'),'SimOptions'))
        prjpath = fileparts(getSimOption('prjname'));
        if ~isempty(prjpath)
            paths{end+1} = prjpath;
        end
    end      
    paths = unique(paths,'stable');
    [j,k] = meshgrid(1:numel(pattern),1:numel(paths));
    pattern = cellfun(@(p,f) fullfile(p,f),paths(k),pattern(j),'unif',0);
    
    files = unique(pickfile(pattern,Inf,'-fullpath'),'stable');    
end

function [files_A,cover_A] = mergedcover(files_A,cover_A,files_B,cover_B)
% [files_C,cover_C] = mergedcover(files_A,cover_A,files_B,cover_B)
% Compose cover indices cover_A, and cover_B (referring to file lists files_A, and files_B, resp.)
% into a unified files_C list and cover_C index.

    f = cover_B == 0;
    cover_B(~f) = cover_B(~f) + numel(files_A);
    f = cover_A == 0;
    cover_A(f) = cover_B(f);
    [files_A,ic] = unique([files_A,files_B],'stable');
    f = cover_A > 0;
    cover_A(f) = ic(cover_A(f));
end

function [gaps,msg] = listgaps(missing,UTCinterval)
% Given a boolean vector MISSING corresponding to days UTCINTERVAL(1):UTCINTERVAL(2), generate a
% cell array of couples {[a1,b1],[a2,b2],..} with the start and end times [ai,bi] of each data gap 
% (consecutive zeros in MISSING). Return also a message MSG summarizing the gaps.

    if ~any(missing)
        gaps = {}; msg = 'No missing data'; return;
    end
    
    t = UTCinterval(1):UTCinterval(2)-0.5;
    assert(numel(missing) == numel(t),'Expecting %d element boolean array',numel(t));

    if all(missing)
        gaps = {[t(1),t(end)+1]}; 
        msg = sprintf('Complete data period (%s to %s) missing',datestr(t(1)),datestr(t(end))); 
        return;
    end
    
    chg = diff(missing(:));
    starts = find([missing(1);chg] > 0);

    ends = [chg;-missing(end)] < 0;
    gaps = mat2cell([t(starts),t(ends)+1],ones(numel(starts),1),2);

    if nargin > 1
        msg = arrayfun(@(s,e) [datestr(s) ' to ' datestr(e)],t(starts),t(ends),'unif',0);
        msg = shortliststr(msg,'missing data period');
    end        
end

function UTCinterval = useful_UTC_interval(interval,Loc)
% B = USEFUL_UTC_INTERVAL(A,LOC) - Get a span B of whole UTC days that contains all light hours 
%   of the given interval A. That is:
%
%   a) B.TimeZone = 'UTC', and dateshift(B,'start','day') == B
%   b) Either B(1) < A(1) or sunel(t) < 0 for all t in [A(1),B(1)]
%   c) Either B(2) > A(2) or sunel(t) < 0 for all t in [B(2),A(2)]
%
% Where sunel(t) is the true solar elevation at time t for location LOC.

    UTCinterval = interval(:);
    UTCinterval.TimeZone = 'UTC';
    
    UTCinterval(1) = dateshift(UTCinterval(1),'start','day');
    UTCinterval(2) = dateshift(UTCinterval(2)-eps(1),'end','day');
    
    if diff(UTCinterval) == days(1), return; end
    
    UTCinterval(2) = UTCinterval(2) - 1; % start of last day
    
    YMD = [year(UTCinterval),month(UTCinterval),day(UTCinterval)];
    [~,~,sr,ss] = solarposition.daytime(Loc,YMD);
    
    UTCinterval(2) = UTCinterval(2) + 1; % end of last day
    
    if interval(1) > ss(1), UTCinterval(1) = UTCinterval(1) + 1; end
    if interval(2) < sr(2), UTCinterval(2) = UTCinterval(2) - 1; end
end

function [files,cover] = filtercandidates(files,type,Loc,UTCinterval,dt,opt)
% [F,COVER] = FILTERCANDIDATES(FILES,TYPE,LOC,UTCINTERVAL,DT,OPT)
%   Go through a list of FILES, check if they are appropiate candidates for file TYPE, LOC,
%   UTCINTERVAL, and timestep DT, and return a filtered list F of those that seem useful, along
%   with an index COVER, with one element for each day in UTCINTERVAL, and 0 < COVER(j) = k for
%   any UTC day j whose data is to be found in file k.
%
% FILTERCANDIDATES doesn't read all the files. First (unless OPT.ignorenames = true) it attempts 
%   to get details (coordinates and time step) from the file names, assuming they have been 
%   generated by METEOFILENAME (e.g. 'mcclear_21N076_104W608_15m_2018.csv').
%   Files for locations further than OPT.MaxDist, or with time-steps larger (or not integer 
%   divisors) of DT are filtered out quietly.
%
%   On a second pass, it attempts to read the header* of each shortlisted file, using read_CAMS
%   or read_MERRA2 with the flag 'data',false.
%   Files that result on error are dropped quietly if OPT.justguessing, or with a warning
%   otherwise. At this point, if RUNNINGFROMUI, FILTERCANDIDATES will offer to rename any valid 
%   files that don't have a standardized METEOFILE name.
%   Again, files for locations further than OPT.MaxDist, or with time-steps larger (or not integer 
%   divisors) of DT are filtered out quietly.
%
%   Finally, an index COVER is created based on the data interval found on each header. Files with
%   summarization ~ DT, and with the largest overlap are given priority.

    % interval must span an integer number of days in UTC
    t0 = UTCinterval(1);
    n = days(UTCinterval(2)-t0);
    cover = zeros(n,1);
    
    if isempty(files), return; end
    files = files(:)';

    % Search for std. names e.g. 'mcclear_21N076_104W608_PT15M_2018.csv'
    pattern = ['.*' type '_(?<latitude>[\dNS]+)_(?<longitude>[\dEW]+)_(?<dt>P[^_]+)_.*'];
    F = regexpi(files,pattern,'names'); 
    stdz = ~cellfun(@isempty,F); % standardized

    if any(stdz)
        dd2deg = @(s) str2double(regexprep(s,'[NSEW]','.')).*(1-2*contains(s,{'W','S'}));
        for j = find(stdz)
            try
                F{j}.latitude = dd2deg(F{j}.latitude);
                F{j}.longitude = dd2deg(F{j}.longitude);
                F{j}.dt = isoduration(F{j}.dt);
            catch
                stdz(j) = false;
            end
        end
    end

    if ~opt.ignorenames && any(stdz)
    % Filter out files that, just from the title, we know will be useless
        F = structarray(F(stdz),1);
        d = solarposition.arcdist(Loc.latitude,Loc.longitude,[F.latitude],[F.longitude],6.4e6);
        tocheck = ~stdz;
        tocheck(stdz) = d <= opt.MaxDist & mod(dt,[F.dt]) == 0;
        files = files(tocheck);
        stdz = stdz(tocheck);
        if isempty(files), return; end
    end
    
    % Read the header of files that have made the short list
    switch type
    case 'merra2', readfcn = @read_MERRA2;
    case 'mcclear', readfcn = @read_CAMS;
    end
    F = cell(size(files));
    valid = true(size(files));
    for j = 1:numel(files)
        try F{j} = readfcn(files{j},'data',false);
        catch
            valid(j) = false;
        end
    end
    if ~all(valid)
        if ~opt.justguessing
           warning('getremote:badfiles',['Reading the following file(s) resulted on error:\n',...
                    shortliststr(files(~valid),'-newlines')]);
        end
        files = files(valid);
        if isempty(files), return; end
        F = F(valid);
        stdz = stdz(valid);
    end
    F = structarray(F,1);
        
    if any(~stdz) && runningfromUI()
    % offer to rename to standardized names, for faster evaluation next time
        msg = sprintf('Do you want to rename %d %s files for faster filtering next time?',...
            nnz(~stdz),type);
        switch optquestdlg(msg,'getremote:rename','yes','no','no')
        case 'yes'
            newnames = arrayfun(@(F) meteofilename(F,F.date_begin,F.date_end,F.dt,...
                'prefix',type,'suffix','.csv'),F(~stdz),'unif',0)';
            newnames = cellfun(@(old,new) fullfile(fileparts(old),new),files(~stdz),newnames,'unif',0);
            try cellfun(@movefile,files(~stdz),newnames); 
            catch ERR, warning(getReport(ERR)); 
            end
        end
    end
    
    % filter again, this time checking for interval overlaps
    d = solarposition.arcdist(Loc.latitude,Loc.longitude,[F.latitude],[F.longitude],6.4e6);
    tocheck = d <= opt.MaxDist & mod(dt,[F.dt]) == 0 & ...
        UTCinterval(1) < [F.date_end] & UTCinterval(2) > [F.date_begin];
    
    files = files(tocheck);
    F = F(tocheck);
    if isempty(files), return; end
    
    if numel(files) > 1
        overlap = min([[F.date_end] - UTCinterval(1); UTCinterval(2) - [F.date_begin];...
                       [F.date_end] - [F.date_begin]],[],1);
        overlap = min(overlap,UTCinterval(2) - UTCinterval(1));

        % Use first files with dt ~ F.dt, and with the largest overlap
        [~,idx] = sortrows([dt./[F.dt]; days(overlap)]','descend');
        files = files(idx);
        F = F(idx);
    end

    t = t0 +(0:n-1)';
    missing = true(n,1);
    for j = 1:numel(files)
        fromj = false(n,1);
        fromj(missing) = F(j).date_begin <= t(missing) & F(j).date_end > t(missing);
        cover(fromj) = j;
        missing(fromj) = false;
        if ~any(missing)
            files(j+1:end) = [];
            break; 
        end
    end
end

function [files,H] = remote_requests(type,Loc,gaps,dt,opt,varargin)
% Handle WPS requests, split interval (if too long for a single request), and cast errors into
% warnings (lets GETREMOTE only crash at the end, if the ammount of missing data turns out to
% be critical).
    
    if isempty(gaps), files = {}; return; end
    
    optargs = [{'MaxDist',opt.MaxDist,'verbose',opt.verbose},varargin{:}];

    switch type
    case 'merra2'
        
        span = cellfun(@(c) days(c(2)-c(1)),gaps);
        if dt < minutes(10)
        % 1 min and 5 min time-steps are limited to one month of data per delivery
            maxspan = 31;
        elseif dt < days(1)
        % Time-steps below 1 day are limited to one year of data per delivery
            maxspan = 366;
        end
        toolarge = find(span > maxspan);
        for j = toolarge
        % Divide the requested period in retrievable chunks
            n = round(span/ceil(span/maxspan));
            t = gaps{j}(1);
            while t < gaps{j}(2)
                gaps{end+1} = [t,min(gaps{j}(2),t+days(n))]; %#ok<AGROW>
                t = gaps{end}(2);
            end
        end
        gaps(toolarge) =[];
        assert(numel(gaps) <= opt.MaxRequests,'Max. number of partial-requests exceeded.');
    
        reqfcn = @(g) wps_MERRA2(Loc,g,dt,optargs{:});
    case 'mcclear'
        reqfcn = @(g) wps_CAMS_McClear(Loc,g,dt,optargs{:});
    end
    
    n = numel(gaps);
    H = cell(n,1);
    files = cell(n,1);
    working = true(n,1);
    msg = {};
    for j = 1:n
        try
            [H{j},files{j}] = reqfcn(gaps{j});
        catch ERR
            msg{end+1} = ERR.message; %#ok<AGROW>
            working(j) = false;
        end
    end
    if ~all(working)
        msg = uniquecell(msg);
        warning(strjoin(msg,newline()));
        files = files(working);
        H = H(working);
    end
end

function [D,meta] = mergedata(type,files,cover,UTCinterval,dt)
% [D,M] = MERGEDATA(TYPE,FILES,COVER,UTCINTERVAL,DT) - read FILES and merge their data into a 
%   single structure, selecting the data source for each day UTCINTERVAL(1):UTCINTERVAL(2) 
%   according to index COVER, and downsampling data (as required) to timestep DT.
    
    switch type
    case 'merra2'
        D = arrayfun(@read_MERRA2,files,'unif',0);
    case 'mcclear'
        D = arrayfun(@read_CAMS,files,'unif',0);
    end
    N = days(diff(UTCinterval));
        
    M = cell(numel(files),1);
    for j = 1:numel(files)
    % Make sure time-stamps are regular, unique, centered intervals
    
        try
            D{j}.t = parsetime(D{j}.t,'-regular','latticeunit',hours(1),'interval','c','TimeZone','UTC');
        catch ERR
            if debugging(), keyboard();
            else, rethrow(ERR);
            end
        end
        [D{j},M{j}] = struct2timetable(D{j},'t');

        dayidx = floor(days(D{j}.t - UTCinterval(1))) + 1;
        used = dayidx > 0 & dayidx <= N;
        used(used) = cover(dayidx(used)) == j;
        D{j} = D{j}(used,:);
                
        if M{j}.dt ~= dt
            D{j}.t = D{j}.t - M{j}.dt/2; % BOI
            D{j} = retime(D{j},'regular','mean','timestep',dt); % (§)
            D{j}.t = D{j}.t + dt/2; % COI
        end
        
        % NOTE(§): MATLAB's RETIME works here only because we are confident that BOI-summarized
        % time steps of both MERRA2 and MCCLEAR data have labels that match the hour-aligned grid
        % of RETIME (0:00,DT,2DT..1:00,..).
        
        D{j}.file(:) = j; % keep track of data origin
    end

    D = cat(1,D{:});
    D = sortrows(D,'t');
    % D = timetable2struct(D);
    
    if numel(M) == 1
        meta = M{1};
        meta.files = files{1};
        D.file = [];
    else
        M = structarray(M);

        % Copy common metadata directly, keep inconsistent fields as array D.files
        meta.files(numel(M)) = struct();
        [meta.files.filename] = deal(files{:});

        for f = fieldnames(M)'
            x = uniquecell({M.(f{1})});
            if isscalar(x)
                meta.(f{1}) = x{1};
            else
                [meta.files.(f{1})] = deal(M.(f{1}));
            end
        end
    end
end

function [MD,msg] = altitude_correction(MD,h_data,h_site)
% Apply corrections due to altitude differences between MERRA2-pixel-mean and site

    LAPSE_RATE = -0.0065; % K/m
    MIN_OFFSET = 10;      % m
        
    dh = h_site - h_data;
    if all(abs(dh) <= MIN_OFFSET)
        msg = {sprintf('No altitude correction, |offset| <= %0.1f m',max(abs(dh)))};
        return; 
    end

    if isstruct(MD)
        % parsestruct(MD,{'Ta','Patm'},'opt', {'air_density','RH'},'-n','-r','-e');
        % isfield = @isfield;
    elseif istimetable(MD)
        % parsestruct(timetable2struct(MD),{'Ta','Patm'},'opt', {'air_density','RH'},'-n','-r');
        isfield = @(T,x) any(strcmp(T.Properties.VariableNames,x));
    end
    compatiblesize(dh,MD.Ta);
    
    if ~all(isfinite(dh))
        warning('Unknown altitude offset, cannot apply corrections');
        return;
    end
    
    if isfield(MD,'RH'), RH = MD.RH; else, RH = 0.5; end
    assert(~any(RH < 0 | RH > 1.2),'Unexpected humidity, wrong units?');
    
    if isfield(MD,'air_density'), rho = MD.air_density;
    else, rho = wetair.density(MD.Ta,RH,MD.Patm);
    end
    assert(~any(rho < 0.4 | rho > 2),'Unexpected density, wrong units?');

    Rspec = MD.Patm./(rho.*(MD.Ta+273.15));                  % specific gas constant, wet air
    Ta_fix = MD.Ta + LAPSE_RATE*(dh);
    P_fix = MD.Patm.*((MD.Ta+273.15)./(Ta_fix+273.15)).^(9.806./(LAPSE_RATE*Rspec));

    msg = {sprintf('Applying %+0.0f m altitude correction (site: %0.0f, original data: %0.0f):',...
            dh,h_site,h_data)};
        
    trackchange('Temperature',MD.Ta,Ta_fix,'°C');
    trackchange('Pressure',MD.Patm/100,P_fix/100,'hPa');

    if isfield(MD,'air_density')
        rho_fix = rho.*((MD.Ta+273.15)./(Ta_fix+273.15)).^(9.806./(LAPSE_RATE*Rspec));
        trackchange('Density',rho*1000,rho_fix*1000,'g/m³');
        MD.air_density = rho_fix; 
    end
    
    if isfield(MD,'RH')
        w = wetair.mixratio(wetair.vpres(MD.Ta,RH),MD.Patm);     % mix-ratio [kg/kg]
        RH_fix = wetair.relhum(Ta_fix, w, P_fix); 
        trackchange('Humidity',MD.RH*100,RH_fix*100,'%');
        MD.RH = RH_fix;
    end

    MD.Ta = Ta_fix; 
    MD.Patm = P_fix; 

    function trackchange(var,xo,xn,unit)
        msg{end+1,1} = sprintf('%s: %+0.2f %s (std. %0.2f %s) on mean %0.1f %s',....
                                var,mean(xn-xo),unit,std(xn-xo),unit,mean(xo),unit);
    end
end

function D = sync(D,t,dt,idx)
% D = SYNC(D,T,DT,IDX) - Resample, downsample, and/or offset data D to match timesteps T. 
%   D is expected to be a timetable with regular timesteps, labels at interval centers.
%   T is either a regular time vector with step DT (labels at center of interval), or a time  
%   vector of instant samples (sorted and unique, not necessarily regular), signaled by DT = 0.
%   A sorting/repetition index IDX can be used to return D.x(IDX,:) for each variable D.x

    if nargin < 4, idx = ':'; end
    % assert(istimetable(D) && ~isempty(D.Properties.TimeStep));

    if dt > 0
        offset = mod(t(1) - dt/2 - dateshift(t(1),'start','hour'),dt)/dt;
        if offset > 0.5, offset = offset - 1; end

        if D.Properties.TimeStep ~= dt || offset ~= 0
           m = minutes([D.Properties.TimeStep,dt]);
           D = resamplestructure(D,m,'-centered','offset',offset);
           
           % numerical precission can account for errors of order << milliseconds
           D.t = datetime(round(datevec(D.t)*1000)/1000,'TimeZone',D.t.TimeZone);
        end
        D = retime(D,t,'fillwithmissing');
    else
        D = retime(D,t,'linear');
    end

    D = D(idx,:);
    
    % Remove variables with no data (quietly)
    for j = 1:size(D,2)
        if isnumeric(D.(j)) && ~any(isfinite(D.(j))), D(:,j) = []; end
    end
end

