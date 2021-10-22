function [F,data,flaglist,msg] = minQC(data,varargin)
% F =  MINQC(X, 'Name', Value,...) - Detect and flag invalid data, outliers, "flats", and 
%   optionally, observations with low availability or close to other flagged points. Data set X 
%   can be an [N,M] matrix, table, or timetable, representing N observations on M variables.
%
%   F will be an N-vector of integers, encoding 7 binary flags:
%
%       dec2bin(FLAG) = [ system, buffer, MAH, IQR, flat, num, NA ]
%
%   1 NA - missing-value function [ default @isnan ] retured 0.
%   2 num - invalid-value function [ default @(x) isfinite(x) & isreal(x) ] returned 1. 
%
%   4 flat - local variability is under expected threshold: s = MOVMAD(x,W,1)/0.6745 < S, for
%           'window' size W (default 5) and a threshold S defined by 'flatstd' and 'flatmethod'.
%
%   8 IQR - By default, flags values beyond the Inter-Quartile-Range estimated bounds for each
%           variable: Q1 - K·(Q3-Q1) < x < Q3 + K·(Q3-Q1), where Q1 and Q3 are the first and third
%           quartiles for each column, and the factor K is set by the outlier probability 'P',  
%           [from GETSIMOPTION('outliers.P')], namely: K = {norminv(1-P/2)/0.6745-1}/2.
%
%           More generally, a custom quantile 'Q' (e.g. 0.1 for interdecile range) can be set 
%           instead of the default Q = 0.25 (interquartile range). The  scale factor K then becomes
%           K = (norminv(1-P/2)/norminv(1-Q)-1)/2. The '-asym' flag allows assymetric estimates 
%           2·(Q-mu) on each side of the median mu. 'minrange' can be used to guarantee a known
%           minimum spread on each side of the median.
%
%           Alternatively, if the variable(s) follow a known, strongly non-normal distribution D,
%           ('fitdist') parametric confidence bounds can be estimated using FITDIST(X,D).
%
%   [16 MAH] - Mahalanobis distance of observation (based on ROBUSTCOV estimate) is beyond outlier
%           probability: (X-mu)/C·(X-mu) > chi2inv(1-P,M), AND the value seems (at least partly)
%           responsible. Namely, one of the following is true:
%             A. Removing the value results in (X-mu)/c·(X-mu) < chi2inv(1-R,M-1), where R >= P
%                (default R = P) is a predefined confidence bound for "reasonable" observations. 
%             B. The blame (as in A) cannot be assigned to other variable(s) in the observation.
%           The frequency of B cases (complete flagged observations) increases with R, while 
%           reducing the chances of false negatives (Type II errors).
%
%           NOTE: this test is not very robust for non-gaussian distributions, so it is not
%           performed by default. Use 'sequence','all' or 'sequence',{..,'MAH',..} to include.
%
%   32 buffer - point is within 'buffer' steps of a point flagged for other reasons. 
%           Default is M (deactivated).
%
%   64 system - there are more than 'maxinvalid' variables flagged for the same observation.
%           Default is Inf (deactivated).
%
%   Simultaneous flag combinations are encoded into numbers 0-127 (see FLAGLIST output below), 
%   e.g. F(j,k) = 13 = 1 + 4 + 8 implies that variable k at observation j was parsed invalid, 
%   and flagged by both IQR and MAH.
%
%   Specific Quality-Checks and their behavior can be tuned by {Name, Value} pairs (shown below 
%   with their default values):
%
%     'sequence',C - By default, tests will be run in the order set above. The sequence can
%           be defined arbitrarily by a vector of test keys {'num','flat','IQR',..}, the
%           keyword 'all' or list-complements, e.g. {'^','MAH'} (default), as recognized by
%           PARSELIST. Steps can be omitted or even repeated, if desired.
%
%     'independent',false - By default, values flagged by one test are set to NaN before execution
%           of the next. Setting 'independent',true overrides this behavior, making all tests use
%           all data points (sequence becomes irrelevant).
%
%     'valid', @(x) isfinite(x) & isreal(x) - Function handle with [N,M] boolean output.
%   
%     'window', 5 - moving half-window size for flat detection using MOVMAD
%     'flatstd', 0 - hard threshold for flat detection. Interpretation given by 'flatmethod'
%     'flatmethod','IQR-log' - free string containing keys {abs,rel,loc,glo,iqr,log}.
%
%         Strings like: 'local-relative' set a hard threshold for s = MOVMAD(x,W,1)/0.6745 < S
%         based on 'flatstd',K. Namely:
%
%           'abs' xor 'rel': defines whether S* is absolute (S = K) or relative (S = KS*)
%           'loc' xor 'glo': for 'rel', defines whether S* is the global or moving-window STD
%        
%         Strings containing 'IQR' (and optionally 'log') additionally pass MOVMAD estimates for
%         local STD to an assymetric IQR outlier detection algorithm, effectively setting a non-
%         parametric threshold S. Since the distribution for s (chi-square-ish?) is heavily
%         skewed, the IQR filter works best when applied on log(s) instead of s.
% 
%     'P', SimOptions.outliers.P - probability threshold for outliers (typ. 1e-6)
%     'Q',0.25 - Use a custom quantile Q for inter-quantile-range confidence bounds.
%     'asym',false - Use assymetric estimates 2·(Q3-mu) and 2·(mu-Q1) instead of (Q3-Q1).
%     'minrange',0 - Use min(S,Q3-Q1) to set a minimu spread for the data.
%      
%     'R', 0.1 - confidence bound for "reasonable" observations (MAH)
%     'mu', median(X,1) - overrides estimate of central moment. It can be used to force IQR and
%           MAH to work with strongly non-normal distributions (e.g. one-sided). 
%     'cov', ROBUSTCOVARIANCE(X) or MAD(X,1)/0.6745 if this fails - overrides covariance estimate,
%           used both for MAH and for flat detection, as S = sqrt(diag(C))
%
%     'IQRk', {norminv(1-P/2)/norminv(1-Q)-1}/2  - override IQR scale factor
%     'chi2q', chi2inv(1-[P,R],M) - override chi² scores for MAH
% 
%     'maxinvalid', M - max .flagged variables in a row, before flagging the whole row
%     'buffer', 0 - flag any values within 'buffer' observations of already flagged values
%
%     
    % FLAGS might contain keys that are unknown to MINQC (e.g. from METEOQC), however, all keys in
    % OPT.sequence should be in MINQC_FLAGS.
%
% [F,Xr,FLAGS,MSG] = MINQC(X,..) - Return also a copy Xr of X in which flagged points have been 
%       removed (replaced by NaNs); the list of flag keys FLAGS (ordered LSB-MSB), and a
%       summary MSG of the results, of the type: 
%
%           'X/MN values in Y/N observations flagged by QC (Z% missing, W% invalid, ... )'
%
% See also: BINNEDIQR, COMPLETEMETEODATA, CHECKUTCOFFSET

    if nargin == 0, minQC_test(); return; end
    
    MINQC_FLAGS = {'NA','num','flat','IQR','MAH','buffer','system'}; 

    % Get [N,M] matrix from whatever is in DATA
    if isa(data,'timetable')
        X = table2array(data);
    elseif isa(data,'table')
        X = table2array(data);
    elseif isnumeric(data) && ismatrix(data)
        X = data;
    end
    [N,M] = size(X);
    
    [OPT,varargin] = getflagoptions(varargin,{'-independent','-asym','-plot'});

    OPT.sequence = {'^','MAH'};
    OPT.valid = @(x) isfinite(x) & isreal(x);
    OPT.missing = @isnan;
    OPT.mu = [];
    OPT.cov = [];
    OPT.P = getSimOption('outliers.P');
    OPT.Q = 0.25;
    OPT.R = [];
    OPT.IQRk = []; 
    OPT.minrange = 0;
    OPT.fitdist = '';
    OPT.chi2q = [];
    %OPT.flatstd = 1e-4;
    %OPT.flatmethod = 'local-relative';
    OPT.flatstd = 0;
    OPT.flatmethod = 'IQR-log';
    OPT.window = 5;
    OPT.maxinvalid = M;
    OPT.buffer = 0;
    OPT.flaglist = MINQC_FLAGS; 

    OPT = getpairedoptions(varargin,OPT,'restchk');
    
    parsestruct(OPT,{'sequence','flaglist'},'class',{'char','cell','string'},'nonempty');
    [OPT.sequence,seqidx] = parselist(OPT.sequence,MINQC_FLAGS,'flag-key','in test sequence');
    [~,bitidx] = parselist(MINQC_FLAGS,OPT.flaglist,'-soft');    
    
    seqbits = bitidx(seqidx);
    if ~all(seqbits > 0)
       undefined = uniquecell(OPT.sequence(seqbits == 0));
       error('No key in flag list for %s', shortliststr(undefined,'test')); 
    end
    flaglist = OPT.flaglist;
    flaglist(bitidx(bitidx > 0)) = MINQC_FLAGS(bitidx > 0);

    if isempty(OPT.IQRk)
        OPT.IQRk = 0.5*(norminv(1-0.5*OPT.P)/norminv(1-OPT.Q)-1);
    end
    if isempty(OPT.chi2q)
        OPT.chi2q = chi2inv(1-OPT.P,M);
    end
    if isempty(OPT.R), OPT.R = OPT.P; end
    
    parsestruct(OPT,{'window','buffer','maxinvalid'},'-n','-r','-s','-i','-p');
    parsestruct(OPT,{'flatstd','IQRk','chi2q','minrange'},'-n','-r','-s','-p');
    parsestruct(OPT,{'P','Q','R'},'-n','-r','-s','>=',0,'<=',1);
    OPT.chi2q(2) = chi2inv(1-OPT.R,M-1);     

    F = zeros([N,M],meteoQC.TYPE);

    if OPT.independent && any(contains(OPT.sequence,{'flat','IQR','MAH'}))
        [mu,C,s] = updatemoments(X,OPT);
    end
    
    for k = 1:numel(OPT.sequence)
        
        if ~OPT.independent && contains(OPT.sequence{k},{'flat','IQR','MAH'})
            [mu,C,s] = updatemoments(X,OPT);
        end
            
        switch OPT.sequence{k}
        case 'NA'
            bad = OPT.missing(X);
        case 'num'
            bad = ~OPT.valid(X);
        case 'flat' 
            bad = detectflats(X,OPT,s);
        case 'IQR'
            if ~isfinite(OPT.IQRk), continue; end
            bad = IQRoutliers(X,mu,OPT);
        case 'MAH'
            bad = MCDoutliers(X,C,mu,OPT.chi2q);
        case 'buffer'
            if OPT.buffer == 0, continue; end
            bad = buffer(F > 0,OPT.buffer);
            % bad = bad & (F == 0);
        case 'system' 
            if OPT.maxinvalid >= M, continue; end
            bad = sum(F > 0,2) > OPT.maxinvalid;
            bad = repmat(bad,1,M);
        otherwise
            error('You really should not be here');
        end
        
        F(bad) = bitset(F(bad),seqbits(k));
        
        if ~OPT.independent, X(bad) = NaN; end
    end
    
    if OPT.plot, QCplot(data,F,mu,s,flaglist); end

    if nargout > 1
        X(F > 0) = NaN;
        if isa(data,'timetable') || isa(data,'table')
            data{:,:} = X;
        else
            data = X;
        end
    end
    
    if nargout > 2, msg = flagsummary(F,MINQC_FLAGS); end
end

function [mu,C,s] = updatemoments(X,OPT)
% Estimate mean and covariance, or use overrides if included in OPT
            
    M = size(X,2);
    mu = mean(X,1,'omitnan');
    s = std(X,1,1,'omitnan');
    constant = s < sqrt(eps(mu));
    C = zeros(M);
    if all(constant), return; end

    if isempty(OPT.cov)
        full = all(isfinite(X),2);
        if nnz(full) < 2*M
            % robustcov needs sample size >= twice the number of variables.
            warning('minQC:cov',['Not enough full observations to estimate covariance ',...
                                'using uncorrelated variable MAD']);
            s(~constant) = mean(abs(X-mu),1,'omitnan')/0.6744897501960818;
            C = diag(s.^2);
        else
            [C(~constant,~constant),mu(~constant)] = robustcov(X(full,~constant));
            s = sqrt(diag(C))';
        end
    else
        C = OPT.cov;
        if isscalar(C), C = C*eye(M); end
        assert(isequal(size(C),[M,M]),'Expecting [M,M] covariance override');
        s = sqrt(diag(C))';
    end
    if ~isempty(OPT.mu)
        mu = OPT.mu;
        if isscalar(mu), mu = repmat(mu,1,M); end
        assert(isequal(size(mu),[1,M]),'Expecting [1,M] mean override');
    end
        
end

function flagged = detectflats(X,OPT,S)
% Flag 'flats', i.e. periods with abnormally low MOVMAD, 

    KEYS = {'abs','rel','loc','glo','iqr','log'};

    OPT.flatmethod = lower(OPT.flatmethod);
    MTD = cell2struct(cellfun(@(x) contains(OPT.flatmethod,x),KEYS','unif',0),KEYS);
    
    if OPT.flatstd == 0, MTD.abs = 1; MTD.rel = 0;
    else
        assert(xor(MTD.abs,MTD.rel) & xor(MTD.loc,MTD.glo),'Incompatible options');
    end
    
    s = movstd(X,OPT.window,0,1,'omitnan','endpoints','fill');
    % s = movmad(X,OPT.window,1,'omitnan','endpoints','fill')/0.6744897501960818;
    
    if MTD.iqr || (MTD.rel && MTD.loc)
    % Get an estimate of the expected MOVSTD in samples with good availability
        a = movsum(isfinite(X)*1,OPT.window,1,'omitnan','endpoints','fill');
        r = s;
        r(a < OPT.window/2) = NaN;
        mu = median(r,1,'omitnan');
    end
    
    if MTD.abs
        S = OPT.flatstd;
    else
        if MTD.loc, S = mu; end
        S = S.*OPT.flatstd;
    end
    flagged = s <= S;
    
    if MTD.iqr
        OPT.asym = true;
        OPT.fitdist = '';
        if MTD.log
            [lo,~] = IQRoutliers(log(s),log(mu),OPT);
        else
            [lo,~] = IQRoutliers(s,mu,OPT);
        end
        flagged = flagged | lo;
    end
    
    % flagged = buffer(flagged,window);
end

function flagged = buffer(flagged,b)
% Flag b right- and -left neighbors of any already flagged point
    
    if b > 0
        xt = [repmat(flagged(1,:),b,1);flagged;repmat(flagged(end,:),b,1)];
        for j = [-b:-1,1:b]
           xt = xt | circshift(xt,j,1); 
        end
        flagged = xt(b+1:end-b,:);
    end
end

function varargout = IQRoutliers(X,mu,OPT)
% Flag outliers based on IQR (symmetric or assymetric), or on known distribution

    if ~isempty(OPT.fitdist)
    % Fit distribution opt.fitdist to each data column, and calculate confidence bounds
    
        [N,M] = size(X);
        pd = fitdist(X(:),OPT.fitdist,'by',repelem((1:M)',N));
        lo = cellfun(@(pd) icdf(pd,OPT.P),pd,'unif',0);
        hi = cellfun(@(pd) icdf(pd,1-OPT.P),pd,'unif',0);
        lo = [lo{:}];
        hi = [hi{:}];
    else  
    % Inter-Quantile-Range estimate
    
        K = OPT.IQRk;
        q = quantile(X,[OPT.Q,1-OPT.Q],1);

        % Estimate data limits
        if OPT.asym
            lo = min(mu,q(1,:)) - 2*K*max(OPT.minrange/2,mu - q(1,:)); % Q1 - 2 K(mu-Q1)
            hi = max(mu,q(2,:)) + 2*K*max(OPT.minrange/2,q(2,:) - mu); % Q3 + 2 K(Q3-mu)
        else
            IQR = K*max(OPT.minrange,q(2,:)-q(1,:)); % [Q1 to Q3]    
            lo = min(mu,q(1,:)) - IQR;
            hi = max(mu,q(2,:)) + IQR;
        end
    end
    if nargout == 2, varargout = {X < lo,X > hi};
    else, varargout{1} = X > hi | X < lo;
    end
end

function [flagged,mah] = MCDoutliers(X,C,mu,chi2quantile)
% Flag outliers based on Mahalanobis distance

    [N,M] = size(X);
    
    % Fill gaps with medians, and recalculate Mahalanobis distance, to allow flagging of
    % incomplete observations.
    X = X - mu;
    X(~isfinite(X)) = 0;
    mah2 = dot(X/C,X,2);
    out = mah2 > chi2quantile(1);
    
    if M < 2, flagged = out; return; end
        
    % If removing variable j of an outlier observation brings the observation back to 
    % normal, point the blame at j
    flagged = false(N,M);
    for j = 1:M
        f = (1:M) ~= j;
        % x = X(out,:);
        % x(:,j) = x(:,f)*(C(f,f)\C(f,j)); % conditional E(xj) given all other x
        % mhj2 = dot(x/sig,x,2);
        x = X(out,f);
        mhj2 = dot(x/C(f,f),x,2);
        flagged(out,j) = mhj2 < chi2quantile(2);
    end
    
    % If nobody confessed, everyone is guilty!
    allguilty = out;
    allguilty(out) = ~any(flagged(out,:),2);
    flagged(allguilty,:) = true;
    
    if nargout > 1, mah = sqrt(mah2); end
end

function msg = flagsummary(F,FLAGS)
% Generate a summary string of the type: 
%   'X/MN values in Y/N observations flagged by QC (Z% missing, W% invalid, ... )'

    [N,M] = size(F);
    
    Nf = nnz(F > 0);
    if Nf == 0, msg = 'No values flagged by QC'; return; end
    
    msg = sprintf('%d/%d values ',Nf,N*M);
    if size(F,2) > 1
        r = any(F > 0);
        msg = sprintf('%sin %d/%d observations ',msg,nnz(r),N);
    end
    msg = [msg 'flagged by QC'];
   
    % Expand F bit-wise, so that F(i,j,k) is the jth bit of observation i, variable k
    % B = fliplr(dec2bin(F(:),numel(FLAGS)) == '1'); % LSB to MSB
    B = bitget(repmat(F(:),1,numel(FLAGS)),repmat(1:numel(FLAGS),N*M,1));
    B = permute(reshape(B,N,M,[]),[1,3,2]);

    share = sum(B,[1,3])/Nf;
    if nnz(share) == 1
        msg = sprintf('%s (%s)',msg,FLAGS{share > 0});
    else
        share = arrayfun(@(j) sprintf('%0.1f%% %s',share(j)*100,FLAGS{j}),find(share > 0),'unif',0);
        msg = sprintf('%s (%s)',msg,strjoin(share,', '));
    end
end

function QCplot(data,F,mu,s,FLAGS)
% Plot series of standardized variables, with markers for flagged observations

    MINQC_FLAGS = {'NA','num','flat','IQR','MAH','buffer','system'}; 
    MARKERS = {'','x','s','v','^','d','o'};
    
    [~,b] = parselist(MINQC_FLAGS,FLAGS);
    if ~any(b), return; end
    
    FLAGS = MINQC_FLAGS(b > 0);
    b = b(b > 0);
    MARKERS = MARKERS(b > 0);

    if isa(data,'timetable')
        t = data.Properties.RowTimes;
        X = table2array(data); 
        datatags = cellfun(@(v) arrayfun(@(k) sprintf('%s_%d',v,k),1:size(data.(v),2),'unif',0),...
                data.Properties.VariableNames,'unif',0);
        datatags = cat(2,datatags{:});
    elseif isa(data,'table')
        X = table2array(data);
        t = (1:size(X,1))';
        datatags = data.Properties.VariableNames;
    elseif isnumeric(data) && ismatrix(data)
        X = data;
        t = (1:size(X,1))';
        datatags = arrayfun(@(k) sprintf('col. %d',k),1:size(X,2),'unif',0);
    end
    [~,M] = size(X);
    
    T = repmat(t,1,M);
    X = X - mu;
    Z = X./s;
    % Z = sqrt(dot(X/C,X.2));
    
    GUIfigure('minQC','minQC','3:1'); clf(); hold on;
    h_data = plot(T,Z);

    h_flags = gobjects(numel(FLAGS),1);
    msg = cell(1,numel(FLAGS));
    used = false(numel(FLAGS),1);
    for j = 1:numel(FLAGS)
        bad = bitget(F,b(j)) > 0;
        if ~any(bad), continue; end
        used(j) =  true;
        h_flags(j) = plot(T(bad),Z(bad),MARKERS{j});
        msg{j} = sprintf('%s: %d/%d points (%0.2f%%)\n',...
                            FLAGS{j}, nnz(bad),numel(bad),nnz(bad)/numel(bad)*100);
    end
    h_flags = h_flags(used);
    msg = msg(used);
    
    % Remove flagged points from lines
    Z(F > 0) = NaN;
    for j = 1:M
        set(h_data(j),'YData',Z(:,j));
    end
    
    % Don't clutter legend with more than 5 data-tags
    if M > 5
        h_data(5:end) = [];
        datatags(5:end) = [];
        datatags{4} = '...';
    end

    legend([h_flags;h_data],[msg,datatags]);
end

function minQC_test()
% Performance test, for debugging and default-parameter tuning, uses a VARM model with random
% parameters, and random ammounts of noise ("invalid" points, flats, and outliers) - see RUNTEST
% below for details. Performs the analysis RUNS times, plots a histogram of Type I and Type II
% detection errors, and re-runs worst cases with -plot flag for display.

    M = 3;      % variables
    N = 1000;   % observations (per run)
    RUNS = 50;
    
    FLAGS = {'foo','NA','num','flat','bar','IQR','MAH','buffer','system','bam'}; 
    
    SEQ = 'all';

    OPTS = {'P',0.0001,'sequence',SEQ,'-independent','flaglist',FLAGS};
    
    [~,~,FLAGS] = minQC(rand(3,1),OPTS{:});
    
    TypeI = false([N,M]);
    TypeII = false([N,M]);
    tampered = false([N,M]);
    F = zeros([N,M]);
    
%     tampered_mah = false([N,M]);

    rng(1234);
    RNG = randi(100,RUNS,1);
    err = zeros(RUNS,2);
    for k = 1:RUNS
        rng(RNG(k));
        runtest();
        err(k,1) = nnz(TypeI)./nnz(~tampered);
        err(k,2) = nnz(TypeII)./nnz(tampered);
    end
    
    GUIfigure('minQC_test'); clf();
    corrplot(err/nnz(tampered),'Varnames',{'TypI','TypII'});
    
    [~,interesting(1)] = max(err(:,1)); % worst Type I
    [~,interesting(2)] = max(err(:,2)); % worst Type II
    [~,idx] = sort(sum(err,2),'descend');
    interesting(3:5) = [1,floor(RUNS/2),RUNS];   % worst, avg, best
    [interesting,ic] = unique(interesting,'stable');
    
    tags = {'worst_I','worst_II','worst','avg','best'};
    tags = tags(ic);
    
    for k = 1:numel(interesting)
        rng(RNG(idx(interesting(k))));
        runtest('plot',true);
        addtitle(tags{k});
    end

    function lst = flaglist(X) 
        X = unique(X(:)); 
        n = numel(FLAGS);
        B = bitget(repmat(X,1,n),repmat(1:n,numel(X),1)) > 0;
        lst = arrayfun(@(j) strjoin(FLAGS(B(j,:)),' & '),1:size(X,1),'unif',0);
        lst = shortliststr(lst,'flag');
    end
    
    function addtitle(tag)
        
        ax = findobj(get(GUIfigure('minQC'),'children'),'type','Axes');
        fig = GUIfigure(['minQC_test_' tag]);
        nuax = copyobj(ax,fig);
                
        if any(TypeI)
            msg{1} = sprintf('%d false outliers (%0.1f%%), with %s',...
                             nnz(TypeI),nnz(TypeI)/(N*M)*100,flaglist(F(TypeI)));
        else
            msg{1} = sprintf('No false outliers');
        end
        if any(TypeII)
            msg{2} = sprintf('%d/%d uncaught outliers',nnz(TypeII),nnz(tampered));
        else
            msg{2} = sprintf('No uncaught outliers');
        end

        % % Mahalanobis check
        % b = find(contains(FLAGS,'MAH'));
        % MAH = bitget(F,b) > 0;
        % TypeI_MAH = MAH & ~tampered_mah;
        % TypeII_MAH = ~MAH & tampered_mah;
        % if any(TypeI_MAH)
        %     msg{3} = sprintf('%d MAH false outlier observations (%0.1f%%)',nnz(TypeI_MAH),nnz(TypeI_MAH)/(N*M)*100);
        % else
        %     msg{3} = sprintf('No MAH false outlier observations');
        % end
        % if any(TypeII_MAH)
        %     msg{4} = sprintf('%d/%d MAH uncaught outlier observations',nnz(TypeII_MAH),nnz(tampered_mah));
        % else
        %     msg{4} = sprintf('No uncaught MAH outlier observations');
        % end
        
        title(nuax,msg);
    end
        
    function runtest(varargin)
        
        c = rand(M,1)*3;
        AR = {rand(M)-0.5};
        S = cov(rand(M+1,M));
        Mdl = varm('Constant',c,'AR',AR,'Covariance',S);
        X = simulate(Mdl,N);

        tampered = false([N,M]);
        
        % fault states = {missing, invalid, outliers, flats}
        a = ([1 1 1 1] + [10,10,10,10].*rand(1,4))/100;  % asymptotic frequency
        t = [1 1 1 5] + [10,5,1,20].*rand(1,4);  % mean fault duration (geometrical distribution)
        
        % Simulate faults as a Markov-process (independent channels, exclusive faults) 
        R = 1./t;               % probability of repair (jump to state 0)
        f = a.*R/(1-sum(a));    % probability of failure (jump to state j)
        P = [1-sum(f),f;R',diag(1-R)];
        S = simulate(dtmc(P),N-1,'x0',M*eye(1,numel(t)+1)) - 1; % S = 0 = working
        
        % missing value
        bad = S == 1;
        % bad = randperm(numel(X),randi(N/2));
        X(bad) = NaN;
        tampered(bad) = true;

        % "invalid" value
        bad = S == 2;
        % bad = randperm(numel(X),randi(N/2));
        X(bad) = 0.999;
        tampered(bad) = true;

        % far outliers
        bad = S == 3;
        % bad = randperm(numel(X),randi(N/4));
        X(bad) = (X(bad)+1)*2;  
        tampered(bad) = true;
%         tampered_mah = bad;
        
%         T = repmat((1:N)',1,M);
%         GUIfigure('debug'); clf(); hold on;
%         plot(T,X);
%         plot(T(bad),X(bad),'ro');

        % flats (all channels)
        bad = S == 4;
        m = find(diff(bad) > 0);
        if ~isempty(m)     
            X(bad) = interp1([m;M*N+1],[X(m);0],find(bad),'previous');
            tampered(bad) = true;
        end
        % bad = randi(N/2) + (1:randi(N/10));
        % m = randi(M);
        % X(bad,m) = X(bad(1),m);
        % tampered(bad,m) = true;

        F = minQC(X,'valid',@(x) x ~= 0.999,OPTS{:},varargin{:});
        TypeI = ~tampered & (F > 0);  % false positives
        TypeII = tampered & (F == 0); % false negatives
    end
end
