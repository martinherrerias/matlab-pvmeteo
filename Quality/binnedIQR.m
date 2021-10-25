function varargout = binnedIQR(x,y,varargin)
% [XB,LO,HI,E] = BINNEDIQR(X,Y,[N]) - X-binned Inter-Quantile-Range estimate of Y-data spread.
%   Classifies the data into N-1 bins with equal number of points, and estimates the limits of
%   'reasonably expected' data points at each edge XB(j) based on the Inter-Quantile-Ranges of
%   the Y values in the neighboring bins. That is:
%   
%       LO(j) = Q1(j) - K·(QN-Q1)
%       HI(j) = QN(j) + K·(QN-Q1)
%       E = Y > HI | Y < LO
%       where Q1,QN are the Q and (1-Q) quantiles of Yj = Y(XB(j+1) < X < XB(j+1))
%       and K is such that, for a normal distribution P(LO) = P(HI) = P
%
% BINNEDIQR(..,'P',P) - Use a custom probability P, instead of the default SimOptions.outliers.P.
% BINNEDIQR(..,'Q',Q) - Use a custom quantile Q (e.g. 0.1 for interdecile range) instead of the
%   default Q = 0.25 (interquartile range).
% BINNEDIQR(..,'-asym') - Use assymetric estimates 2·(QN-M) and 2·(M-Q1) instead of (QN-Q1).
% BINNEDIQR(..,'minrange',S) - Use min(S,QN-Q1) to set a minimu spread for the data.
% BINNEDIQR(..,'xlim',[A,B]) - Extrapolate (nearest neighbor) LO and HI to ensure that XB covers
%   the interval [A,B].
%
% [L,H,E] = BINNEDIQR(..,'-interp') - Return interpolants L(x) and H(x) that can be evaluated to
%   estimate the lower and upper limits for (new) arbitrary X.
%
% [P,E] = BINNEDIQR(..,'-poly') - Return a closed polygonal envelope instead of the raw points.  
%   The polygon can be used directly to test for outliers, e.g. E = ~INSIDEPOLYGON(P,X,Y). This
%   flag sets a non-zero minrange = 1e-6 as default, to avoid degenerate polygons. 
%
% P = BINNEDIQR(..,'-poly','ends',KEY) - The ends of the polygonal envelope are by default
%   closed at the maximum and minimum X (or at hard limits set by XLIM). They can instead be 
%   rounded (KEY = 'round') or extended (KEY = 'square') with radius/offset (HI-LO)/2.
%
% P = BINNEDIQR(..,'fitdist',DIST) - [Experimental] for heavily skewed-distributed data, quantile
%   estimates for dispersion are not very reliable (even with the -asym option). If the type of
%   expected distribution for each bin is known, FITDIST(Yj,DIST) is applied to every bin, and the
%   extremes LO, HI calculated analytically from the fitted distribution.
%
%   TODO: replace by ISOUTLIER with 'movmethod'??
%
% See also: IQR, FITDIST, INSIDEPOLYGON, KTKDFILTER, ISOUTLIER

    [opt,varargin] = getflagoptions(varargin,{'-poly','-asym','-unif','-interp','-plot'});
    
    opt.N = min(10,ceil(numel(x)/20));
    opt.ends = 'butt';
    opt.minrange = opt.poly.*1e-6;
    opt.xlim = [];
    opt.nrends = 8;
    opt.fitdist = '';
    opt.Q = 0.25;
    opt.P = getSimOption('outliers.P');
    
    [opt,varargin] = getpairedoptions(varargin,opt);
    if ~isempty(varargin)
        assert(numel(varargin) == 1,'Unrecognized argument(s)');
        opt.N = varargin{1};
    end
    if isempty(opt.P)
       if ~isempty(opt.fitdist), opt.P = 0.01; else, opt.P = 0.25; end 
    end
    parsestruct(opt,{'N','minrange','P'},'-n','-r','-p','-s','-f');
    assert(ischar(opt.fitdist),'Invalid opt.fitdist, expecting distribution name');
    N = opt.N;
    
    assert(nnz([opt.poly,opt.interp]) <= 1,'Inconsistent output flags: -poly / -interp');
    
    if ~isempty(opt.xlim)
        x(x < opt.xlim(1) | x > opt.xlim(2)) = NaN;
    end

    if opt.unif
        xb = linspace(min(x,[],'omitnan'),max(x,[],'omitnan'),N)';
    else
    % bin data (irregular bin widths, ~same number of points per bin)
        xb = quantile(x,linspace(0,1,N))';
        xb = unique(xb);
        N = numel(xb);
    end
    b = discretize(x,xb);

    if ~isempty(opt.fitdist)
    % Fit distribution opt.fitdist to each bin (first and last) or to each bin-pair (left and right
    % of each xb) for intermediate breaks. This is done by counting each point twice, once labeled
    % at b and then at b-1:
        
        pd = fitdist([y;y],opt.fitdist,'by',[b;b-1]);
        yb_lo = cellfun(@(d) icdf(d,opt.P),pd)';
        yb_hi = cellfun(@(d) icdf(d,1-opt.P),pd)';
    else  
    % get median, 1st and 3rd quartiles for each bin / bin pair
    
        % IQR-ranges from Q,1-Q to reach P at each side of a normal distribution
        K = (norminv(1-opt.P)/norminv(1-opt.Q)-1)/2; 
    
        q = zeros(N,3);
        for j = 1:N
            yj = y(b == j | b == j-1);
            q(j,:) = quantile(yj,[opt.Q 0.5 1-opt.Q]);
        end

        % Estimate data limits
        if opt.asym
            s = 2*K*max(opt.minrange/2,diff(q,1,2)); % [Q1 to median, median to Q3]    
        else
            s = K*max(opt.minrange,q(:,3)-q(:,1)); % [Q1 to Q3]    
            s(:,2) = s;
        end
        yb_lo = q(:,1) - s(:,1); 
        yb_hi = q(:,3) + s(:,2);
    end
    
    % if ~isempty(opt.fitdist)
    % % Map back to original distribution
    %     f = @(y) icdf(opt.fitdist,normcdf(y),pd{:});
    %     yb_lo = f(yb_lo);
    %     yb_hi = f(yb_hi);
    % end
    
    % Extend to xlimits, if required
    if ~isempty(opt.xlim) && xb(1) > opt.xlim(1)
        xb = [opt.xlim(1);xb];
        yb_lo = [yb_lo(1);yb_lo];
        yb_hi = [yb_hi(1);yb_hi];
        N = N+1;
    end
    if ~isempty(opt.xlim) && xb(N) < opt.xlim(2)
        xb = [xb;opt.xlim(2)];
        yb_lo = [yb_lo;yb_lo(N)];
        yb_hi = [yb_hi;yb_hi(N)];
        N = N+1;
    end
    out = ~inpolygon(x,y,[xb;flipud(xb)],[yb_lo;flipud(yb_hi)]);

    if opt.poly
    % create outline (note first and last points extended up to 0 and 1.2 kt)
    
        r = (yb_hi - yb_lo)/2;
        switch lower(opt.ends)
            case 'butt'
            case 'square'
                xb = [xb(1)-r(1);xb;xb(N)+r(N)];
                yb_lo = [yb_lo(1);yb_lo;yb_lo(N)];
                yb_hi = [yb_hi(1);yb_hi;yb_hi(N)];
            case 'round'
                th = (0.5:opt.nrends)'*90/opt.nrends;
                xb = [xb(1)-r(1)*cosd(th);
                      xb;
                      xb(N)+r(N)*sind(th)];
                yb_lo = [yb_lo(1)+r(1)*(1-sind(th));
                         yb_lo;
                         yb_lo(N)+r(N)*(1-cosd(th))];
                yb_hi = [yb_hi(1)-r(1)*(1-sind(th));
                         yb_hi;
                         yb_hi(N)-r(N)*(1-cosd(th))];
            otherwise
                error('Unknonwn opt.ends');
        end
        xb = [xb;flipud(xb)];
        yb = [yb_hi;flipud(yb_lo)];

        varargout = {polygon(xb,yb),out};
    elseif opt.interp
    % return interpolants   
        L = griddedInterpolant(xb,yb_lo,'linear','nearest');
        H = griddedInterpolant(xb,yb_hi,'linear','nearest');
        varargout = {L,H,out};
    else
        varargout = {xb,yb_lo,yb_hi,out};
    end
    
   
    if opt.plot
        GUIfigure('binnedIQR'); clf(); hold on;
        densityplot(x,y,2);
        scatter(x(out),y(out),2,'r');
        if opt.poly
            polyplot(varargout{1},[0 1 0 0.2],'g');
        else
            plot(xb,yb_lo,':')
            plot(xb,yb_hi,':');
        end
    end
end