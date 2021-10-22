function plot_metrics(ign,BOT,es,bins,lbl,BOTlbl)

    if nargin < 4 || isempty(bins), bins = 20; end
    if nargin < 5 || isempty(lbl), lbl = ''; end
    if nargin < 6, BOTlbl = 'Box density Ordinate Transform'; end

    GUIfigure('knkd_density_test_metrics','','3:1');  % clf(); 

    subplot(1,3,1); hold on; 
    if isinf(bins)
        pd = fitdist(ign(isfinite(ign)),'kernel','Kernel','normal');
        ii = linspace(-10,10,500);
        plot(ii,pdf(pd,ii));
    else
        histogram(ign,bins,'normalization','pdf','DisplayStyle','stairs');
    end
    xlabel('Ignorance Score');
    ylabel('pdf'); 
    grid on; axis square; % xticks(-10:5:10); xlim([-10,10]);

    subplot(1,3,2); hold on;
    if isinf(bins)
        pd = fitdist(BOT(isfinite(BOT)),'kernel','Kernel','normal');
        ii = linspace(0,1,200);
        plot(ii,pdf(pd,ii));
    else
        histogram(BOT,bins,'normalization','pdf','DisplayStyle','stairs');
    end
    xlabel(BOTlbl);
    ylabel('pdf');
    grid on; axis square; xticks(0:0.25:1); xlim([0,1]);

    subplot(1,3,3); hold on;
    if isinf(bins)
        pd = fitdist(es(isfinite(es)),'kernel','Kernel','normal');
        ii = linspace(0,0.6,200);
        plot(ii,pdf(pd,ii));
    else
        histogram(es,bins,'normalization','pdf','DisplayStyle','stairs');
    end
    xlabel('Energy Score');
    ylabel('pdf');
    grid on; axis square; % xticks(0:0.15:0.6); xlim([0,0.6]);

    if ~isempty(lbl)
        lbl = regexprep(lbl,'\$(.*= )(.*)\$','\$$2\$');
        
        lgd = get(gca,'legend');
        if isempty(lgd), lgd = lbl; else, lgd = [lgd.String(1:end-1),lbl]; end
        legend(lgd,'interpreter','latex','fontsize',12,'box','off');
    end
end