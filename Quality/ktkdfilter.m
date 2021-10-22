function varargout = ktkdfilter(kt,kd,varargin)
% OUT = KTKDFILTER(KT,KD,MODE) - Classify each pair Kt(j), Kd(j) of clearness-index vs diffuse-
%   fraction values as outlier (OUT(j) = true) or not, depending on the criteria set by MODE.
%   By default, the criteria is a mixture of predefined polygonal envelope (samples from [1]) (*)
%   and binned statistical analysis.
%
% [OUT,P] = KTKDFILTER(KT,KD,MODE) - return also the polygonal envelope used as outlier criterion.
%
% .. KTKDFILTER(..,'-silent') - Unless '-silent' flag is used, KTKDFILTER calls KTKDPLOT,
%   highlights the outliers, and draws the polygonal envelope of normality.
%
% .. KTKDFILTER(..,'offset',S) - apply an offset of S to predefined envelope. OFFSET is also used
%       as the minimum range (q3 -q1) to ensure envelope from quartile analysis is not degenerate.
% .. KTKDFILTER(..,'nbins',N) - Use N bins for quantile analysis of outliers.
% .. KTKDFILTER(..,'Q',Q) - Quantile to use for BINNEDIQR (see 'sta' method below), Def. Q = 0.02.
% .. KTKDFILTER(..,NAME,VAL) - Provide additional arguments to BINNEDIQR.
%
% INPUT:
%   KT, KD: equal-sized numerical arrays of clearness-index and diffuse-fraction values.
%   MODE: character key (only first three characters are checked):
%     'pre' - use predefined polygonal envelope for normality. Based on the two examples from [1],
%        hand-fitted envelopes for [2], and an optional OFFSET (see below).
%     'sta' - binned analysis of outliers based on quantiles using BINNEDIQR(kn,Y,..). Note that
%        filtering is performed on a transformed space: kn = kt·(1-kd), Y = kt·(1+3·kd)/2. The
%        quantity kt·(1+3·kd)/2 has no physical meaning, but seems to be near-orthogonal to
%        kn.
%     'mix' - (Default) use a mixture of both methods above (a point is qualified as an outlier 
%       only if it falls outside both the predefined and the quantile envelopes).
%
% OUTPUT: 
%   OUT - boolean array of size(KT), where OUT(j) = true means point (j) is not normal.
%   P - POLYGON object, envelope used as outlier criterion.
%
% REFERENCES:
% [1] T. Muneer and F. Fairooz, “Quality control of solar radiation and sunshine measurements –
%   lessons learnt from processing worldwide databases,” Building Services Engineering Research 
%   and Technology, vol. 23, no. 3, pp. 151–166, Aug. 2002.
%
% [2] Marion, B., Anderberg, A., Deline, C., Cueto, J. del, Muller, M., Perrin, G., Rodriguez, 
%   J., Rummel, S., Silverman, T.J., Vignola, F., Kessler, R., Peterson, J., Barkaszi, S., Jacobs, 
%   M., Riedel, N., Pratt, L., King, B., 2014. New data set for validating PV module performance 
%   models, in: 2014 IEEE 40th Photovoltaic Specialist Conference (PVSC). pp. 1362–1366.
%
% FUTURE: review more robust/recent/extended practices on:
% [3] S. Younes, R. Claywell, and T. Muneer, “Quality control of solar radiation data: 
%   Present status and proposed new approaches,” Energy, vol. 30, no. 9, pp. 1533–1549.
%
% See also: KTKDPLOT, COMPLETEMETEODATA

    [opt,varargin] = getflagoptions(varargin,{'-silent'});
    
    opt.nbins = 16;
    opt.Q = 0.02;
    opt.offset = 0;
    
    [opt,varargin] = getpairedoptions(varargin,opt);
    
    if ~isempty(varargin) && ischar(varargin{1})
        mode = varargin{1}; 
        varargin = varargin(2:end);
    else
        mode = 'auto';
    end
    
    sz = size(kd);
    assert(isnumeric(kt) && isnumeric(kd) && isequal(size(kd),sz),...
        'Expecting equal-sized numerical arrays Kd, Kt');
    kt = kt(:); kd = kd(:);
    
    switch lower(mode(1:3))
        case {'aut','mix','hyb'} % statistical outliers & outside predefined envelope
            p = mergepolygons(predefined_env(opt),IQR_envelope(kt,kd,opt,varargin{:}));
            
        case {'pre','env'}, p = predefined_env(opt);      % predefined envelope

        case {'sta','qua'}, p = IQR_envelope(kt,kd,opt,varargin{:});  % statistical outliers
            
        case 'man', error('Interactive mode is no longer supported')
        otherwise, error('Unrecognized operation mode');
    end
    
    outliers = ~insidepolygon(p,kt,kd);

    if ~opt.silent
        ktkdplot(kt,kd,outliers);
        polyplot(p,'none','c','linestyle',':');
    end
    
    if nargout > 0, varargout{1} = reshape(outliers,sz); end
    if nargout > 1, varargout{2} = p; end
end

function p = abs_env()
% Absolute limits: (0 < GHI < 1.2 GNE) & (0 < DHI < 0.8 GNE) recommended by CIE
% equivalent to (0 < kt < 1.2) & (0 < kdkt < 0.8)
    p = polygon([1.2,0,0,0.8:0.05:1.2;
                   0,0,1,0.8./(0.8:0.05:1.2)]');
end

function p = predefined_env(opt)
% Predefined 'normality' envelope (ideally based on a large data set of good data)

    % Known (site specific) envelopes:
    
    % [1] Samples from  T. Muneer and F. Fairooz
    XY{1} = [0,1;0.49,1;1,0.24;1,0;0.89,0;0.11,0.89;0,0.89;0,1];  % Bahrain
    XY{2} = [0,1;0.5,1;1,0.14;1,0;0.8,0;0.2,0.9;0,0.9;0,1];       % UK
    
    % [2] NREL's data for Validating Models (estimated to fit 99.8% of points, hand simplified)
    XY{3} = [0,1;0.57,1;0.98,0.53;1.1,0.51;1.1,0.42;0.8,0.077;0.73,0.085;0.42,0.32;0.18,0.61;0.1,0.92;0,0.97]; % Cocoa
    XY{4} = [0,1;0.57,1;0.83,0.8;1,0.47;1,0.21;0.84,0.057;0.77,0.062;0.34,0.37]; % Golden
    XY{5} = [0,1;0.64,1;0.89,0.69;1,0.47;1,0.26;0.81,0.052;0.77,0.053;0.45,0.25;0,0.93]; % Eugene
    % ... add further (reliable) envelopes here
  
    p = cellfun(@polygon,XY,'unif',0);
    p = mergepolygons([p{:}],'pos','pos'); 
    p = offsetpolygon(p,opt.offset);
    p = intersectpolygons(p,abs_env());
end

function p = IQR_envelope(kt,kd,opt,varargin)
% Statistical normality envelope, based on binned outlier analysis

    e = -0.1; % doesn't help for filtering, but required for back-transformation @ kt = 0

    X = kt.*(1-kd); % kn 
    % Y = kt.*(1+kd);
    Y = kt.*(1+3*kd)/2;  % no physical meaning, but approx. orthogonal to kn
    
    p = binnedIQR(X,Y,opt.nbins,'-poly','-asym','Q',opt.Q,varargin{:});
    %c = abs_env();
    %c = polygon(c.x.*(1-c.y),c.x.*(1+3*c.y)/2);
    %p = intersectpolygons(p,c);
    % p = refinepoly(p,0.05);
    
    % kt_b = (p.x + p.y)/2;        % kt·(1-kd) + kt·(1+kd) = 2·kt
    % kd_b = (p.y-p.x)./(2*kt_b);  % kt·(1-kd) - kt·(1+kd) = 2·kt·kd
    
    b = e+2*p.y+3*p.x;
    c = 2*p.y-p.x;
    kd_b = (b - sqrt(max(0,b.^2-4.*e.*c)))/(2*e);
    kd_b(c == 0) = b(c == 0)/e;
    kt_b = (b - e*(1+kd_b))/4;
    
    p = intersectpolygons(polygon(kt_b,kd_b),abs_env());
end

% % TODO: Change to a maximum-likelihood approach, that not only discards outliers using sharp
% %   thresholds, but detects (and fixes) precision errors:
% %
% %  a) modify ktkdplot to return kernel-density interpolant D = ktkdplot(kt,kd,bad,args...); 
% %     and/or estimate predefined kernel-density interpolant D based on high-quality data.
%
%   ALPHA = 0.05;
%   N = 100;
%   
%   bad = D(kt,kd) < ALPHA;
% 
% if any(bad,'all')
%     
%    % b) introduce (correlated) noise into kt,kd points. User must provide U_GHI = UNC.GHI./GHI and
%    %    U_DHI = UNC.DHI./DHI
%    noise = randn(2,N);
%    kt_n = kt(bad).*noise(1,:).*U_GHI(bad);
%    kd_n = kd(bad).*(noise(2,:).*U_DHI(bad) - noise(1,:).*U_GHI(bad));
% 
%    % Get the joint likelihood of P & D, D might beed some normalization.
%    P = normcdf(-abs(noise))*2;
%    P = P(1,:).*P(2,:);
%    P = D(kt_n,kd_n).*P; 
% 
%    % c) Points with max joint likelihood P < ALPHA are actually outliers, in  other cases, tweak
%    %    kd,kt input to the values that maximize P & D
%    [P,i] = max(P,[],2);
%    ok = P > ALPHA;
%    i(ok) = sub2ind(size(kd_n),find(ok),i(ok)); % abs. index in kd_n, kt_n
% 
%    tweaked = bad;
%    tweaked(bad) = ok;
%    kd(tweaked) = kd_n(i);
%    kt(tweaked) = kt_n(i);
% 
%    bad = bad & ~tweaked;
% end
    
