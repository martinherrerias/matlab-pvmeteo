function [error,P0,Pb,msg] = checkUTCoffset(GHI,t,varargin)
% OFFSET = CHECKUTCOFFSET(GHI,T,LOC) 
% OFFSET = CHECKUTCOFFSET(GHI,T,CSGHI) 
%   Perform a cross-correlation analysis between GHI and the estimated/provided clear-sky-GHI 
%   according to time-stamps T (datenum, UTC), to detect a possible offset in the time-stamps 
%   (Time-Zone and/or interval-summarization). Offset detection is accomplished in several steps:
%
%       [0. Calculation of CSGHI, using PVLMOD_CLEARSKY_INEICHEN(LOC,T)]
%       1. (quiet) hourly-offset detection (*) & correction, using full data set
%       2. Clear-sky-sample detection using provisional offset
%       3. (Interactive) offset detection (*), using clear samples
%       4. (optional) half time-step offset detection (*), using clear samples
%
%   (*) Offset-detection comprises cross-correlation of GHI vs CSGHI to find the max. correlation
%   lag -- more strictly speaking, the maximum of a lower-bound estimate for correlation,
%   considering the changing statistical sample size -- see FINDLAG for details.
%
% CHECKUTCOFFSET(..,'-fractional') - By default, cross-correlation is performed only with lags
%   of -12:+14 hours ± 1 interval, i.e. assuming that errors are due to either wrong (integer)
%   time-zone, or wrong interval summarization. '-fractional' flag allows searching along the
%   complete interval -12·n:14·n, where n is the number of samples per hour.
%
% CHECKUTCOFFSET(..,'halfstep',true/false) - as a final step, resample data at half-intervals, 
%   to check for wrong interval summarization. Default is true if dt > 5 min.
%
% CHECKUTCOFFSET(..,'-quiet') - By default, whenever a non-zero OFFSET is detected, warnings are
%   issued, and if the session runs interactively, the user is asked to confirm whether to use or
%   discard the offset. '-quiet' mode overrides this behavior, issuing no warnings, and choosing 
%   automatically to keep an offset when P(offset ~= 0) > MINP.
%   
% The following name,value pairs control de details of the algorithm:
%
% ..,'CSGHI',[] - override clear-sky model: provide a size(GHI) time series for comparison
%
% ..,'clearsky',[] - provide boolean vector identifying clear-sky steps, otherwise estimated 
%               from PVLMOD_DETECT_CLEAR_TIMES.
% ..,'minCShours',2 - min. clear-hours per day (otherwise full day is discarded)
% ..,'minCSdays',10 - min. equivalent days (accumulated clear-sky periods, excluding days 
%   discarded due to 'minCShours') required to use clear-sky filter. If PVLMOD_DETECT_CLEAR_TIMES 
%   fails to find enough clear sky points, CHECKUTCOFFSET will issue a warning and use all points.
%
% ..,'minGHI',20 -  threshold values to filter input
% ..,'maxlag' - hours(12)
%    
% ..,'minP',0.95 - minimum P value to consider a correlation improvement significant
%   
% ..,'minoffset',minutes(6) - don't try to detect smaller offsets
% ..,'minstep',minutes(1) - downsample time series above this resolution 
%
% INPUT:
%   GHI: vector of global-horizontal-irradiance values, can contain NaNs
%   t: datenum vector of UNIFORM time-steps, UTC
%   Loc: location structure, as required by PVL_SPA
%
% OUTPUT: DURATION scalar. OFFSET > 0 means GHI is ahead of CSGHI(t), i.e. t_fix = t - OFFSET;
%
% TODO: Seasonal offset detection (insert breaks due to DST!) 
%
% See also: FINDLAG, COMPLETEMETEODATA, PVLMOD_CLEARSKY_INEICHEN

    if nargin == 0, test(); return; end % DEBUG
    narginchk(3,Inf);
    
    assert(isnumeric(GHI) && isvector(GHI) && isreal(GHI),'Invalid GHI');
    assert(numel(t) == numel(GHI),'Inconsistent GHI and t');
    assert(numel(t) > 4,'Time series too short for any assessment')
    
    [opt,varargin] = getflagoptions(varargin,{'-quiet','-halfstep','-fractional'});
    if ~opt.halfstep, opt.halfstep = []; end
    opt.minCSdays = 10;  % equivalent clear-days
    opt.minCShours = 2;  % clear-hours per day (to keep the day)
    opt.minGHI = 20;     % threshold values to filter input
    opt.maxlag = 12; 
    opt.minP = 0.95;     % minimum P value to consider a correlation improvement significant
    opt.minoffset = 6;   % minutes, don't try to detect smaller offsets
    opt.clearsky = [];
    opt.interval = 'c';
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    assert(isscalar(varargin),'Unrecognized syntax');
    if isstruct(varargin{1})
        Loc = parselocation(varargin{1});
        CSGHI = []; % clear-sky model override
    else
        CSGHI = varargin{1};
        validateattributes(CSGHI,{'numeric'},{'real','vector','size',size(GHI)},'','CSGHI'); 
    end

    if ~isempty(opt.clearsky)
        validateattributes(opt.clearsky,{'numeric','logical'},{'binary','vector'});
        if isscalar(opt.clearsky)
            if opt.clearsky, CS = []; else, CS = 0; end
        else
            CS = opt.clearsky > 0;
            validateattributes(CS,{'logical'},{'size',size(GHI)});
        end
    else
       CS = [];
    end

    [t,dt,idx] = parsetime(t,'-regular','interval',{opt.interval,'c'});
    opt.timezone = t.TimeZone; t.TimeZone = 'UTC';
    assert(mod(3600,seconds(dt)) == 0,'Time-step needs to be an integer fraction of 1 hour');
    dt = minutes(dt);

    if isempty(opt.halfstep)
        opt.halfstep = dt >= opt.minoffset*2; 
    end
    validateattributes(opt.maxlag,{'numeric','duration'},{'real','scalar','positive','nonzero'},'','maxoffset');
    validateattributes(opt.minoffset,{'numeric','duration'},{'scalar','integer','positive','nonzero'},'','minoffset');
    if isa(opt.maxlag,'duration'), opt.maxlag = hours(opt.maxlag); end
    if isa(opt.minoffset,'duration'), opt.minoffset = minutes(opt.minoffset); end

    if opt.halfstep && dt/2 < opt.minoffset
       opt.halfstep = false;
       warning('checkutc:minoffset','OPT.minoffset < dt/2, overriding OPT.halfstep');
    end

    % Make sure timesteps are uniform, and comprise an integer number of days
    u = t(end)+minutes(dt):minutes(dt):t(1) + ceil(days(t(end) + minutes(dt) - t(1))) - minutes(dt)/2;
    t = [t;u'];
    if ~isequal(idx',1:numel(t))
        GHI = revertfilter(GHI,idx,size(t),NaN);
        if ~isempty(CSGHI), CSGHI = revertfilter(CSGHI,idx,size(t),NaN); end
        if numel(CS) > 1, CS = revertfilter(CS,idx,size(t),false); end
    end
    
    if dt < 1
    % downsample anything above 1 step/minute (performance)
        [t,GHI,CSGHI] = cellresample(60/dt,60,t,GHI,CSGHI);
        dt = 1;  
    elseif opt.halfstep      
        t = resamplestructure(t,2);
        GHI = interpn(GHI,1)';
        if ~isempty(CSGHI), CSGHI = interpn(CSGHI,1)'; end
        dt = dt/2;
    end

    if isempty(CSGHI)
        K = floor(dt/2)*2;
        if K > 1   
            % Evaluate clear-sky model at ~1 min resolution, to reduce averaging error
            tt = resamplestructure(t,K,'-centered');
            CSGHI = pvlmod_clearsky_ineichen(Loc,tt);
            CSGHI = avgdownsample(CSGHI,K);
        else
            CSGHI = pvlmod_clearsky_ineichen(Loc,t);
        end
    end
    
    GHI(GHI < opt.minGHI) = NaN;
    CSGHI(CSGHI < opt.minGHI) = NaN;
    
    if isempty(CS)
    % (Quietly) detect any major offset, just for clear-sky-detection to work
        if mod(30,dt) == 0
            k = floor(opt.maxlag*2)*30/dt;
            coarse_lags = -k:30/dt:k;
        else
            k = floor(opt.maxlag)*60/dt;
            coarse_lags = -k:60/dt:k;
        end

        [coarse_offset,c,~,lb,P,ub] = findlag(CSGHI,GHI,coarse_lags,'P',opt.minP,'-periodic');

        if P(coarse_lags == 0) < opt.minP, offset = 0; else, offset = coarse_offset; end

        % Select clear samples [using provisional offset correction]
        CS = pvlmod_detect_clear_times(GHI,circshift(CSGHI,offset),dt);
    end
    
    if numel(CS) > 1
        % Selecte clear days (those with at least opt.minCShours)
        [~,~,d] = unique(dateshift(t-minutes(offset*dt),'end','day'));
        idx = accumarray(d,CS) > opt.minCShours*60/dt; 
        CS = idx(d) & CS;

        if nnz(CS)*dt/1440 < opt.minCSdays
            warning('checkutc:nocs','Not enough clear-samples, using all points for UTC-offset check');
            CS = [];
        else
            GHI(~CS) = NaN;
            GHI = GHI(idx(d));
            CSGHI = CSGHI(idx(d)); 
            t = t(idx(d));
        end
    end
    
    if opt.fractional && dt < 60
    % allow anything within range (typ. ± 0:12 hours)
        k = floor(opt.maxlag*60/dt);
        lags = -k:k; 
    else
    % only integer hours in specified range (typ. ± 0:12 hours) ... 
        k = floor(opt.maxlag)*60/dt;
        lags = -k:60/dt:k;
        
        if opt.halfstep, ots = -2:2;            % include -1,-1/2, +1/2 and +1 steps 
        elseif dt > opt.minoffset, ots = -1:1;  % or just ± 1 steps
        else, ots = [];
        end
        if ~isempty(ots)
            lags = reshape(lags + ots',1,[]); 
            lags(abs(lags)*dt/60 > opt.maxlag) = [];
            lags = unique(lags,'sorted');
        end
    end
    
    if ~all(CS) || ~all(ismember([lags,coarse_offset],coarse_lags))
    % Re-evaluate offset-detection for all lags [using clear-sky samples]
    
        [offset,c,~,lb,P,ub] = findlag(CSGHI,GHI,lags,'P',opt.minP,'-periodic');
    else
        lags = coarse_lags;
        offset = coarse_offset;
    end
    P0 = P(lags == 0);
    Pb = P(lags == offset);
        
    error = minutes(offset*dt);
    error.Format = 'hh:mm:ss';
    
    if error == 0
        msg = sprintf('zero offset detected, P(offset~=0) < %0.1f%%',P0*100);
    else
        msg = sprintf('detected offset = %s: P(offset=0) = %0.1f%%',char(error),(1-P0)*100);
        % if ( Pb < opt.minP )
            [~,nextbest] = max(c - 1*(lags == offset));
            nextbest = minutes(lags(nextbest)*dt);
            nextbest.Format = 'hh:mm:ss';
            msg = sprintf('%s\n%0.1f%% confidence vs offset = %s',msg,Pb*100,char(nextbest));
        % end
    end
    
    if error == 0 && P0 >= opt.minP, return; end % confident that error = 0
    
    if ~opt.quiet
        hfig = offset_plot(t,GHI,CSGHI,lags,offset,c,lb,ub,P,dt,msg);   
    end
    
    if error == 0 || opt.quiet || ~runningfromUI()
        if error == 0
            warning('checkutc:unclear',['Unclear ' msg]);
        elseif ( P0 < opt.minP )
            error(:) = 0;
            warning('checkutc:ignore',['Ingoring ' msg]);
        else
            warning('checkutc:apply',['Applying ' msg]);
        end
        return;
    end

    % At this point: error ~= 0 &&  ~opt.quiet && runningfromUI   
    if P < opt.minP, DEF = 'No'; else, DEF = 'Yes'; end

    while ishandle(hfig)
        reply = optquestdlg(sprintf('%s (see plot). Apply correction?',msg),...
            'UTC-check','Yes','No','Wait',DEF);
        if ~strcmp(reply,'Wait'), break; end

        btn = uicontrol('Parent',hfig,'Style','pushbutton','String','Back',...
            'Units','normalized','Position',[0.9 0 0.1 0.1],'Visible','on',...
            'Callback',@(~,~) uiresume());
        uiwait(hfig);
        delete(btn);
    end

    switch reply
        case 'Yes', msg = ['Correcting for ' msg];
        case 'No', msg = ['Discarding ' msg]; error = 0;
        case 'Wait'

        otherwise, error('Aborted by user, wrong UTC offset');
    end
    warning('checkutc:result',msg);
end

function varargout = cellresample(n,m,varargin)
% Resample vectors X, Y, .. (sample rate N) to new sample rate M

    % [varargin{:}] = compatiblesize(varargin{:});
    N = numel(varargin{1});
    
    fld = arrayfun(@(x) num2str(x,'foo_%d'),1:numel(varargin),'unif',0);
    S = cell2struct(varargin(:),fld);
    S = resamplestructure(S,[m,n],N,'-centered'); % [Nm/n,N] resampling matrix
    varargout = struct2cell(S);
    
%     [~,K] = resamplestructure([],[m,n],N,'-centered'); % [Nm/n,N] resampling matrix
% 
%     varargout = varargin;
%     for j = 1:numel(varargin)
%         if isempty(varargin{j}), continue; end % for CSGHI
%         
%         X = double(varargin{j});  % (*)
%         invalid = isnan(X);
%         k = K;
%         k(:,invalid) = 0;
%         X(invalid) = 0;
%         Y = full(K*X)./sum(k,2);
%         varargout{j}(:) = Y(:);   % typecast back (*)
%     end
end

function hfig = offset_plot(t,x,y,lags,offset,c,lb,ub,P,dt,msg)

    hfig = GUIfigure('UTCoffset','UTC-offset check','3:1'); clf(hfig);
    ax = subplot(1,3,1:2); 
    ax2 = subplot(1,3,3); 

    title(ax,msg);
    hold(ax,'on');
    plot(ax,t,x);
    plot(ax,t,y);
    plot(ax,t,circshift(x,-offset));
    ylabel(ax,'GHI, CSGHI (normalized)');
    legend(ax,'GHI','CSGHI','offset GHI','box','off');
    xlim(ax,t(1) + [0,days(2)]);
    plotarrows(ax,'southwest');

    hold(ax2,'on');
    plot(ax2,lags,c,'rx-'); 
    plot(ax2,lags,lb,'r:');
    plot(ax2,lags,ub,'r:');
    
    c(isnan(c)) = 0;
    [~,idx] = sort(c,'descend');
    bestlags = lags(idx(1:4));
    xlim(ax2,[min(bestlags)-1,max(bestlags)+1]);

    xlabel(ax2,sprintf('Lag k (%0.1f min) / unit',dt));
    ylabel(ax2,'Corr( GHI_i, CSGHI_{i+j})');
    grid(ax2,'on');
    
    yyaxis(ax2,'right');
    set(ax2,'yscale','log')
    plot(ax2,lags,P,'o-');
    ylabel('P( corr_{i+k} > corr_{i+j} for best k)');
end

function test()
% TODO: Use realistic synthetic GHI series

%     opt.fractional = false;
%     opt.halfstep = true;
%      
%     opt.minCSdays = 10;  % equivalent clear-days
%     opt.minCShours = 2;  % clear-hours per day (to keep the day)
%     opt.minsamples = 20; % min. overlapping samples to estimate correlation
%      
%     opt.minGHI = 20;  % threshold values to filter input
%     opt.maxlag = [];
%      
%     opt.nDraws = 100; % resampling draws for bootstrapping
%     opt.minP = 0.99;  % minimum P value to consider a correlation improvement significant
%     
%     CSGHI = []; % clear-sky model override
%     
%     opt.minoffset = 3; % minutes, don't try to detect smaller offsets
%     1 = 1;   % minutes, downsample time series above this resolution
    
    [MD,t,Loc] = getMeteoData('*.meteo');
    offset = checkUTCoffset(MD.GHI,t,Loc,'-quiet');
    t = t - offset;
    
    % OFFSET = randi(2)-1;
    OFFSET = randi(5)-4;
    STEPS = randi(2)-1;
    t = circshift(t,STEPS) + hours(OFFSET); % evil laughter
    
    fprintf('Offset = %d hours, %d step\n',OFFSET,STEPS);
    checkUTCoffset(MD.GHI,t,Loc,'minoffset',1,'halfstep',true);
end
