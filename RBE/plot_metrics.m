function plot_metrics(ign,BOT,es,varargin)
    
    opt.bins = 20;
    opt.lbl = '';
    opt.BOTlbl = 'Box Ordinate Transform';
    opt.norm = 'pdf';
    opt.ylim = {[],[],[]}; % IGN, BOT, ES
    opt.xlim = {[],[0 1],[]};
    opt = getpairedoptions(varargin,opt,'dealrest',3);
    
    GUIfigure('knkd_density_test_metrics','','3:1');  % clf(); 
    
    BOT = BOT(isfinite(BOT));
    es = es(isfinite(es));
    ign = ign(isfinite(ign));
    
    X = {ign,BOT,es};
    ax = arrayfun(@(j) subplot(1,3,j),1:3);
    for j = 1:3
        hold(ax(j),'on');
        if isinf(opt.bins)
            pd = fitdist(X{j},'kernel','Kernel','box');
            ii = linspace(pd.icdf(0.001),pd.icdf(0.999),200);
            plot(ax(j),ii,pdf(pd,ii));
        else
            if isscalar(opt.bins) && ~isempty(opt.xlim{j})
               bins = linspace(opt.xlim{j}(1),opt.xlim{j}(2),opt.bins+1);
            else
               bins = opt.bins; 
            end
            histogram(ax(j),X{j},bins,'normalization',opt.norm,'DisplayStyle','stairs');
        end
        switch j
            case 1, xlabel(ax(j),'Ignorance Score');
            case 2, xlabel(ax(j),opt.BOTlbl); xlim([0,1]);
            case 3, xlabel(ax(j),'Energy Score');
        end
        ylabel(ax(j),opt.norm);
        axis(ax(j),'square'); grid(ax(j),'on'); 
        if ~isempty(opt.ylim{j}), ylim(ax(j),opt.ylim{j}); end
        if ~isempty(opt.xlim{j}), xlim(ax(j),opt.xlim{j}); end
    end

    if ~isempty(opt.lbl)
        lbl = regexprep(opt.lbl,'\$(P\(.*\) =)(.*)\$','\$$2\$');
        
        lgd = get(gca,'legend');
        if isempty(lgd), lgd = lbl; else, lgd = [lgd.String(1:end-1),lbl]; end
        legend(lgd,'interpreter','latex','box','off');
    end
end