function TPM = knkd_density_fit(X0,info,TESTS,varargin)
% TPM = KNKD_DENSITY_FIT(X0,INFO,TESTS,..) - Fit a set of Adaptive diffusion Kernel Density
%   Estimates to data X0, as output by KNKD_DENSITY_PREP. 
%
%   TESTS must be a cell-array with elements of one of three kinds:
%
%       TESTS{j} = 'key' - for named keys {'constant','separation','copula','uniform'} creates
%           variations of a constant density estimate (2D interpolants).
%       TESTS{j} = {'a',..} - creates conditional density estimates on variable(s) 'a',.. of X0.
%           These will be 2+N-dimensional interpolants, for N conditioning variables.
%       TESTS{j} = '(a)x(b)x..' - combine pre-calculated models 'a','b',.. These are only
%           'virtual' models, and do not require the creation of a bulky 2+N-D interpolant.
%           See MIXED_MODELS.
%
% Result will be a table TPM with:
%
%    TPM.Properties.RowNames{j} = 'key' / 'a+b' / '(a)x(b+c)'.. name-ID for model/combination as
%           listed above. Note that a+b stands for a model conditioned on both a AND b, whereas
%           (a)x(b) is the uncorrelated product of models a and b.
%
%    TPM.conditions{j} = {}, {'a'}, {'a','b'} ..
%    TPM.interpolant{j} = griddedInterpolant object G, P(kd,kn | x,y,..) = G(kd,kn,x,y..)
%    TPM.texlbl{j} = latex formatted label
%
% Recognized ..,'name',value,.. pairs with their default value:
%
%     'gridtype','a' - 'r' regular/'a' adaptive (quantile-based) evaluation grid
%     'gridsteps',40 - use up to N points for each dimension
%     'uniformshare',0.5 - for adaptive grids, control the mix of quantile-based- and uniform.
%     'minstep',1e-4 - minimum allowed grid gap
%     'gridrange',[0,Inf] - force grids to start/end (non finite values ignored).
%     'gam' = @(N) ceil(sqrt(N)) - function to set AKDE quality (might be reduced due to memory)
%     'plot',true - plot density estimate (or a slice of it)
%     'print',false - save figure copies to disk.
%
% TODO: wrap into a more general class, that allows flexible definition of inputs (currently
%   "conditions") and state variables (currently hardcoded to kn-kd).
%
% See also: KNKD_DENSITY_TEST, KNKD_DENSITY_PREP

    LBL = info.texlbl;
    dt = info.dt;
    N_tests = numel(TESTS);

    VAR = X0.Properties.VariableNames;
    X0 = X0{:,:};
    [N,M] = size(X0);
    assert(isequal(VAR(1:2),{'kn','kd'}),'Expecting kn, kd as first columns');

    opt.gridtype = 'a'; % 'r' regular/'a' adaptive                                  
    opt.gridsteps = 40;
    opt.minstep = 1e-4;
    opt.uniformshare = 0.5;
    opt.gridrange = [0,Inf];
    opt.gam = @(N) ceil(sqrt(N));
    opt.plot = true;
    opt.print = true;
    opt = getpairedoptions(varargin,opt,'restchk');
    
    % LVL = -2:2; % levels to plot
    if opt.print
        if ~isfolder('fig'), mkdir('fig'); end
        printfig = @(name) cellfun(@(fmt) print(gcf,['-d' fmt],'-r800',['./fig/' name '.' fmt]),{'svg','png'});
    else
        printfig = @(name) 0;
    end
    
    if nargin < 3 || isempty(TESTS)
        TESTS = {'uniform','constant',{'cosz'},{'lastkn'},{'lastkd'},{'lastkn','lastkd'}};

        % uncorrelated combinations of above estimates
        MIXES = {'(lastkd)x(lastkn)','(cosz)x(lastkn+lastkd)','(cosz)x(lastkd)x(lastkn)'};
        TESTS = [TESTS,MIXES];
    end
           
    GRIDS = arrayfun(@(j) gridvector(X0(:,j),opt.gridsteps,opt.gridtype,opt.gridrange,...
                                    opt.minstep,opt.uniformshare),1:M,'unif',0);
                  
    g0 = GRIDS(1:2);
    [g0{1},g0{2}] = ndgrid(g0{1:2});
    x0 = [g0{1}(:),g0{2}(:)];

    name = regexprep(info.data,'(^[^_\b]*)(.*)','$1');
    name = sprintf('knkd_density_%s_%s',name,isoduration(dt));
    name = sprintf('%s_%c%du%0.0fgrid',name,opt.gridtype,opt.gridsteps,100*opt.uniformshare);

    if opt.plot
        GUIfigure('knkd_density','',[100,100,500,400]); clf(); hold on;
        plotgrid(g0{1}(:,1),g0{2}(1,:),LBL);
        scatter(x0(:,1),x0(:,2),3,[1,1,1]*0.6,'fill','DisplayName','Grid Points') 
        scatter(X0(:,1),X0(:,2),1,[0.8 0.7 1],'fill','DisplayName','Training data')
        printfig([name '_grid']);
    end

    % Check that all will fit in memmory
    ng = cellfun(@numel,GRIDS);
    d = cellfun(@numel,TESTS)+2;
    Ng = cellfun(@(c) prod(ng(ismember(VAR,[{'kn','kd'},c]))),TESTS);

    % AKDE currently needs max( 2QD, (10+5d)N + 2(N+D²)G ) x 64 bits for Q query points & G = gam.
    assert(maxarraysize() > 2*max(d.*Ng),'Insufficient memory, try reducing grid sizes');
    M = maxarraysize() - max(d.*Ng)/2;
    d = max(d);
    gam = ceil(0.5*(M - (10+5*d)*N)/(N+d^2));
    if isa(opt.gam,'function_handle')
        gam = min(opt.gam(N),gam);
    else
        gam = min(opt.gam,gam);
    end
    fprintf('gam = %d\n',gam);

    name = sprintf('%s_%dgam',name,gam);

    info = ['kn-kd kernel density estimates from meteo data: ',newline(),info.data,...
            sprintf('Adaptive grids at ~%d quantiles, AKDE modes (gam) = %d.\n\n',opt.gridsteps,gam),...
            strjoin([{'Variables:'};info.descriptions],[newline() '  '])];

    TPM.keys = repmat({''},N_tests,1);
    TPM.conditions = TESTS(:);
    TPM.interpolant = cell(N_tests,1);
    TPM.texlbl = cell(N_tests,1);

    % weights for trapezoidal integration
    w = point_weights(GRIDS{1},GRIDS{2});
    
    ismixed = false(N_tests,1);

    for j = 1:N_tests
    try
        if iscell(TESTS{j})
        % Regular case P(x|Y) = f(...)
        
            idx = varindex(VAR,TESTS{j});
            if isempty(TPM.keys{j})
                if numel(idx) <= 2
                    TPM.keys{j} = 'constant';
                    TPM.texlbl{j} = '$P(x_i \mid Y_{i-1}) = P(x_o)$';
                else
                    TPM.keys{j} = strjoin(cellstr(TESTS{j}),'+');
                    TPM.texlbl{j} = ['$P(x_i \mid Y_{i-1}) = f(' strjoin(LBL(idx(3:end)),',') ')$'];
                end
            else
                keyboard();
            end
            fprintf('\nConditions: {%s}\n',TPM.keys{j});

            TPM.interpolant{j} = fitdensity(X0(:,idx),GRIDS(:,idx),w,gam);

        elseif any(strcmp(TESTS{j},{'uniform','copula','constant','separation'}))
        % Named special cases
        
            TPM.conditions{j} = {};
            idx = 1:2;
            
            TPM.keys{j} = TESTS{j};
            switch TPM.keys{j}
            case 'constant'
                TPM.interpolant{j} = fitdensity(X0(:,idx),GRIDS(:,idx),w,gam);
                TPM.texlbl{j} = '$P(x_i \mid Y_{i-1}) = P(x_o)$';
                
            case 'uniform'
                TPM.texlbl{j} = '$P(x_i \mid Y_{i-1}) = 1$';
                TPM.interpolant{j} = griddedInterpolant(g0{1},g0{2},ones(ng(1),ng(2)));
                continue;
                
            case {'copula','separation'}
                
                [~,ie] = ismember('constant',TPM.keys);
                assert(ie > 0 && ~isempty(TPM.interpolant(ie)),...
                    'copula/separation require precalculated "constant" model');
    
                TPM.interpolant{j} = TPM.interpolant{ie};
                P = TPM.interpolant{j}.Values;
                if strcmp(TPM.keys{j},'separation')
                    P = conditionalpdf(P,sum(X0(:,1:2),2),g0{1}+g0{2});
                    TPM.texlbl{j} = '$P(x_i \mid Y_{i-1}) = P(x_o)/P(k_{t,o})$';
                else
                    P = conditionalpdf(P,X0(:,1),g0{1},X0(:,2),g0{2});
                    TPM.texlbl{j} = '$P(x_i \mid Y_{i-1}) = P(x_o)/\prod P(x_{o,j})$';
                end
                P = P./sum(P.*w,1:2);
                TPM.interpolant{j}.Values = P;
            end
        else
        % Just parse to see if TESTS{j} is a valid mix of existing estimates...
            ismixed(j) = true;
            TPM.keys{j} = TESTS{j};
            continue;
        end
        
        if opt.plot
            
            % knkd_density_plot(TPM.interpolant{1}.GridVectors,TPM.interpolant{1}.Values,'pdf',{'k_n','k_d'},'ax',gca(),'ktgrid',[],'sepmdl',{'reindl1'});
            % legend('location','northwest','EdgeColor','none','interpreter','latex')
            printfig([name '_grid']);
        
            knkd_density_plot(TPM.interpolant{j}.GridVectors,TPM.interpolant{j}.Values,'cdf',LBL(idx));
            title(TPM.texlbl{j},'interpreter','latex')
            printfig([name '_' TPM.keys{j}]);
        end

    catch ERR
        TPM.keys{j} = strjoin(cellstr(TESTS{j}),'+');
        warning('TEST: %s failed with error: %s',TPM.keys{j},ERR.message);
    end
        
    end

    TPM = struct2table(TPM);
    TPM.Properties.RowNames = TPM.keys; TPM.keys = [];
    TPM.Properties.Description = info;
    TPM.Properties.UserData.name = name;
    TPM.Properties.UserData.labels = [VAR',LBL];
    
    if any(ismixed)
        TPM(ismixed,:) = [];
        TPM = mixed_models(TPM,TESTS(ismixed));
    end
end

function iv = varindex(VAR,conditions)
    [ic,iv] = ismember([{'kn','kd'},conditions],VAR);
    assert(all(ic),'Unexpected variable names');
end

function G = fitdensity(X,g,w,gam)

    d = size(X,2);
    assert(d == numel(g),'Inconsistent arguments');

    [g{:}] = ndgrid(g{:});
    s = size(g{1});
    
    gg = reshape(cat(d+1,g{:}),[],d);    
    pdfx = akde(X,gg,gam);
    clear gg
    
    pdfx = reshape(pdfx,s);
    cpdfx = pdfx./sum(pdfx.*w,1:2);
    cpdfx(~isfinite(cpdfx)) = 0;

    G = griddedInterpolant(g{:},cpdfx,'linear','nearest');
end

function g = gridvector(x,n,type,range,minstep,k)

    if nargin < 5 || isempty(minstep), minstep = 1e-4; end
    if nargin < 6 || isempty(k), k = 0; end

    if nargin < 4 || isempty(range), range = [Inf,Inf]; end
    [a,b] = bounds(x);
    if isfinite(range(1)), a = range(1); end
    if isfinite(range(2)), b = range(2); end
    
    switch type
    case 'r'
        g = linspace(a,b,n)';
    case 'a'
        s = sort(x);
        u = (s-a)./(b-a);
        p = linspace(0,1,numel(x))'.*(1-k) + k*u;
        if a < s(1), s = [a;s]; p = [0;p]; end
        if b > s(end), s = [s;b]; p = [p;1]; end   
        g = interp1(p,s,linspace(0,1,n)','linear');

        g = uniquetol(g,minstep,'datascale',1);
    end
    
    if isinf(range(1)), g = [2*g(1)-g(2);g]; end
    if isinf(range(2)), g = [g;2*g(end)-g(end-1)]; end
    g = single(g);
end

function cpdfx = conditionalpdf(cpdfx,varargin)

    TINY = 1e-3;

    x = varargin(1:2:end);
    g = varargin(2:2:end);
    
    for j = 1:numel(x)
        pd = fitdist(x{j},'kernel','Kernel','triangle');
        marginal = pdf(pd,g{j});
        cpdfx = cpdfx./(marginal + TINY);
    end
    cpdfx(isnan(cpdfx)) = 0;
    cpdfx = min(cpdfx,prctile(cpdfx(:),99));
end

function plotgrid(x,y,LBL)
    xticks(x); xticklabels(every_nth_lbl(x));
    yticks(y); yticklabels(every_nth_lbl(y));
    xlabel(['$' LBL{1} '$'],'interpreter','latex','fontsize',14);
    ylabel(['$' LBL{2} '$'],'interpreter','latex','fontsize',14);
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    grid on;

    function lbl = every_nth_lbl(x,n)
        if nargin < 2, n = round(numel(x)/4); end
        lbl = repmat({''},1,numel(x));
        vis = mod((1:numel(x))-numel(x),n) == 0;
        lbl(vis) = arrayfun(@(x) num2str(x,'%0.2f'),x(vis),'unif',0);
    end
end
    