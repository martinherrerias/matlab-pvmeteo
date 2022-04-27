function varargout = plot_metrics(X,varargin)
% PLOT_METRICS(X,...) plot histograms of Score metrics X.
% PLOT_METRICS([ign,bot,es,pit.1,pit.2..]) expected order (auto-labeling)
% PLOT_METRICS(X,'typ',{'es','pit','pit',..}) custom columns
%
% ..,'bins',N - bins on each histogram
% ..,'xlim',Q,'ylim',R - Q,R must be a 2 x size(X,2) matrices, with limits for each column
% ..,'lbl',c - set 'DisplayName',c on each plotted line.
% ..,'norm','pdf' - histogram normalization
% ..,'mlbl',C - cell array of custom xlabels for each histogram
% ..,'ax',H - array of axes (e.g. to overlay histograms from other tests)
%
% See also: KNKD_DENSITY_TEST, RBE
   
    if isnumeric(X), X = num2cell(X,1); end
    cellfun(@(x) validateattributes(x,'numeric',{'vector','real'}),X);
    n = size(X,2);
    
    FLD = {'ign','bot','es','pit'};
    MLBL = {'Ignorance Score','BOT','Energy Score','PIT'};
    XLIM = [-Inf,0,-Inf,0; Inf,1,Inf,1];
    YLIM = [-Inf(1,4); Inf(1,4)];

    opt.lbl = '';
    opt.bins = 20;
    opt.typ = FLD; opt.typ(end+1:n) = FLD(end); opt(n+1:end) =[];
    opt.norm = 'pdf';
    opt.ylim = [];
    opt.xlim = [];
    opt.mlbl = {};
    opt.ax = [];
    opt = getpairedoptions(varargin,opt,'dealrest',2);
    validateattributes(opt.lbl,'char',{});
    validateattributes(opt.bins,'numeric',{'integer','positive'});
    
    [opt.typ,idx] = parselist(opt.typ,FLD);
    opt.xlim = parselimits(opt.xlim,XLIM(:,idx));
    opt.ylim = parselimits(opt.ylim,YLIM(:,idx));

    if isempty(opt.mlbl)
        opt.mlbl = uniquelabels(MLBL(idx));
    end
    if isempty(opt.ax)
        r = max(1,floor(sqrt(n)));
        c = ceil(n/r);
        GUIfigure('plot_metrics','',sprintf('%d:%d',c,r)); clf(); 
        ax = arrayfun(@(j) subplot(r,c,j,'tag',opt.mlbl{j}),1:n);
    else
        ax = opt.ax;
        assert(numel(ax) >= n && all(ishandle(ax)),'Bad axes');
    end

    oklimits = @(x) any(isfinite(x)) && ~(x(2) <= x(1));
    for j = 1:n
        hold(ax(j),'on');
        X{j} = X{j}(isfinite(X{j}));
        if ~isempty(X{j})
            if isinf(opt.bins)
                pd = fitdist(X{j},'kernel','Kernel','box');
                if ~isfinite(opt.xlim(1,j)), opt.xlim(1,j) = pd.icdf(0.001); end
                if ~isfinite(opt.xlim(2,j)), opt.xlim(2,j) = pd.icdf(0.999); end
                ii = linspace(opt.xlim(1,j),opt.xlim(2,j),200);
                h = plot(ax(j),ii,pdf(pd,ii));
            else
                if ~isfinite(opt.xlim(1,j)), opt.xlim(1,j) = min(X{j}); end
                if ~isfinite(opt.xlim(2,j)), opt.xlim(2,j) = max(X{j}); end
                if oklimits(opt.xlim(:,j))
                    e = linspace(opt.xlim(1,j),opt.xlim(2,j),opt.bins+1);
                    h = histogram(ax(j),X{j},e,'normalization',opt.norm,'DisplayStyle','stairs');
                else
                    h = histogram(ax(j),X{j},'normalization',opt.norm,'DisplayStyle','stairs');
                end
            end
            if ~isempty(opt.lbl), h.DisplayName = opt.lbl; end
        end
        xlabel(ax(j),opt.mlbl{j});
        ylabel(ax(j),opt.norm);
        axis(ax(j),'square'); grid(ax(j),'on'); 
        if oklimits(opt.ylim(:,j)), ylim(ax(j),opt.ylim(:,j)'); end
        if oklimits(opt.xlim(:,j)), xlim(ax(j),opt.xlim(:,j)'); end
    end

    if nargout > 0, varargout{1} = ax; end
end

function L = parselimits(lim,L)
    if isempty(lim) || ~any(isfinite(L),'all'), return; end
    validateattributes(lim,'numeric',{'size',size(L)});
    set = isfinite(lim);
    L(set) = lim(set);
end

function rlbl = uniquelabels(rlbl)
    [lbl,~,ia] = unique(rlbl);
    rep = accumarray(ia,1) > 1;
    if ~any(rep), return;end
    for j = find(rep)'
        idx = ia == j;
        rlbl(idx) = arrayfun(@(k) sprintf('%s.%d',lbl{j},k),1:nnz(idx),'unif',0);
    end
end