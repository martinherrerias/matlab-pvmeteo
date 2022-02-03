function [Xc,varargout] = fillgaps(X,dt,timescale,varargin)
% Xc = FILLGAPS(X,DT,TSCALE,[VARNAME,MISSING,...])
%   Wrapper for MSSA (Multi-Channel Singular Spectrum Analysis) gap-filling.
%   For improved performance & noise reduction, the signal is first downsampled by TSCALE/DT,
%   gaps filled (MSSA), and the signal again resampled. Non-missing original values are copied
%   replaced in this final series, though the smoothed series Xr can be retrieved as a second 
%   argument (see below).
%
%   X can be a vector or NxD matrix of correlated variables.
%   DT and TSCALE can be DURATION or numeric (hours) scalars.
%
%   ..,'lags' - Default changes with respect to MSSA, to include 1:3-day lags and 1:3-hour lags,
%       on top of the default 1:3-step lags. This is intended to capture daily patterns in the
%       time series.
%
%   ..,'excluded',MISSING - As with MSSA, MISSING values are excluded from gap-filling, which
%       (can improve convergence for large gaps).
%
%   ..,'varname',name - If provided (not empty), a warning will be issued whenever points are 
%   missing.
%
%   ..,'name',val,'-flag',.. Any additional options will be passed directly to MSSA.
%
% [xc,Xr,S,Rc,..] = FILLGAPS(..) - Return original output of MSSA
% 
% See also: MSSA, COMPLETEMETEODATA, BESTESTIMATE

    if nargin == 0, test(); return; end
    narginchk(1,Inf);
    
    opt.varname = '';
    opt.exclude = [];
    opt.lags = [];
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    validateattributes(X,{'numeric'},{'2d','real'},'','X'); 
    N = size(X,1);
    
    if ~isempty(opt.varname)
        validateattributes(opt.varname,{'char'},{},'','varname');
    end

    if isempty(opt.exclude)
        opt.exclude = false(N,1);
    else
        validateattributes(opt.exclude,{'numeric','logical'},{'binary','size',[N,1]},'','exclude');
    end

    available = isfinite(X);
    if all(available | opt.exclude,'all') && nargout < 2
        Xc = X;
        return;
    end
    available = available & ~opt.exclude;
    
    % f = samples per unit-period
    if nargin < 3 || isempty(timescale) || isempty(dt)
        dt = [1,1]*3600;
    else
        if isduration(timescale), timescale = hours(timescale); end
        if isduration(dt), dt = hours(dt); end
        dt = round([dt,timescale]*3600);
    end

    % Downsample if required, according to timescale
    if dt(1) < dt(2)
       % Xr = avgdownsample(X,f,~excluded*1,'full');
       % ex = avgdownsample(excluded*1,f,'full') == 1;
       [Xr,W] = resamplestructure(X,dt,'-centered','-full','-omitnan');
       ex = W*opt.exclude > 1 - 0.5*dt(1)/dt(2);
       % ex = ex & circshift(ex,1) & circshift(ex,-1); % avoid NaNs near gap edges when resampling
    else
       Xr = X;
       ex = opt.exclude;
    end
    
    if isempty(opt.lags)
        % opt.lags = unique(round([1:3/timescale,(1:3)*24/timescale]));
        opt.lags = unique(round([(1:7),(1:7)/timescale,(1:7)*24/timescale]));
        
        [n,D] = size(Xr);
        maxlag = min([n/2,floor(maxarraysize()/(4*n*D))]);
        opt.lags(opt.lags > maxlag) = [];
    end

    % Gap-filling via 
    [Xr,S,RC,L,opt_MSSA] = MSSA(Xr,'exclude',ex,'lags',opt.lags,varargin{:});
        
    % resample, if required
    if dt(1) < dt(2)
        [Xc,W] = resamplestructure(Xr,fliplr(dt),'-centered','-full');
        Xc(N+1:end,:) = [];
        if nargout > 1
            W(N+1:end,:) = [];
            Xr = single(W*double(Xr));
            S = single(W*double(S));
            RC = reshape(RC,size(RC,1),[]);    
            RC = single(W*double(RC));
            RC = reshape(RC,N,size(Xr,2),[]); 
        end
    else
        Xc = Xr;
    end
    Xc(available) = X(available); % Restore original points, whenever valid

    % GUIfigure('debug'); clf(); hold on
    % plot(X); 
    % if f > 1 && nargout <= 1, Xr = W*Xr; end
    % plot(Xr);
    
    if nargout > 1
        varargout = {Xr,S,RC,L,opt_MSSA};
    end
    
    if ~isempty(opt.varname) && ~all(available,'all')
        warning('cmd:missing','interpolating %d/%d missing points in %s',...
            nnz(~available),nnz(~opt.exclude)*size(available,2),opt.varname);
    end
end

function test()

    t = (0.5/1440:1/1440:365)';
    N = numel(t);
    
    % rng(0.123456789);
   
    % some monthly random variation
    m = rand(12,1);
    t0 = solarposition.mdoy(1:12);
    X = interp1([t0(12)-365,t0,t0(1)+365],m([12,1:12,1]),t,'pchip') + 1;

    % ... + an annual cycle
    X = X + (pvl_extraradiation(t+1)/1367-1)/0.05;
    
    % ... + some random harmonics
    L = [1/15,1/7,1/3];
    H = [1/2,1,1.5];
    X = X.*(1 + 0.2*cos(t.*L*pi)*rand(numel(L),1)).*(1 + 0.2*cos(t.*H*pi)*rand(numel(H),1));
    
    % ... and a clear daily pattern
    
    X = X + 0.5*sum([1,0.3*rand(1,2)].*cos(t.*pi).^[1,2,4],2);
    
    X = X - mean(X);
    X = X/std(X);
    Xo = X;
    
    % ... and noise
    NOISE = 0.1;
    noise = randn(N/60,1);
    noise = interp1((0.5/24:1/24:365)',noise,t,'pchip','extrap');
    noise = NOISE*(noise + randn(N,1)*NOISE);
    X = X+noise;

    GAPS = 20;
    GAP_LEN = ceil(0.8*N/GAPS);
    
    d = size(X,2);
    gaps = arrayfun(@(a,n) a + (0:n-1),randi(N*d-GAP_LEN,GAPS,1),randi(GAP_LEN,GAPS,1),'unif',0);
    X([gaps{:}]) = NaN;
    

    [~,Xr,S] = MeteoData.fillgaps(X,minutes(1),hours(1),'-plot');
    GUIfigure('MSSA'); subplot(2,1,2); cla(); hold on;
    plot(Xr+1.96*S,'y-')
    plot(Xr-1.96*S,'y-')
    plot(Xr);
    plot(X,'.');
    plot(Xo);

%     for j = 1:d
%         patch([1:N,N:-1:1]',[Xr(:,j)+1.96*S(:,j);flipud(Xr(:,j)-1.96*S(:,j))],j,...
%             'facealpha',0.2,'edgecolor','none');
%     end

%     Mdl = varm('Constant',NaN(d,1),'lags',[1,1440],'AR',repmat({NaN(d)},1,2));
%     [Mdl,EstSE,logL,E] = estimate(Mdl,X);
    
%     GUIfigure('fillgaps','fillgaps-test','3:1'); clf(); hold on
%     plot(t,X);
%     plot(t,Xc,':');
end
