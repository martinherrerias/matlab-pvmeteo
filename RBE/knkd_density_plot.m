function [h,cb,c] = knkd_density_plot(g,P,varargin)
% KNKD_DENSITY_PLOT(G,P,[SCALE,LBL,AX],...) - plot either 2D contours or 3D iso-surfaces for a 
%   kernel density estimate on kn, kd, [z,..] represented by cell array of grid vectors G and ND 
%   array of PDF values P.
%
%   For 4D+ dimensional estimates, a 3D 'slice' P(:,:,:,j,k,..) with indices j,k,.. near the 
%   median of grids g{4},g{5},..
%
% KNKD_DENSITY_PLOT(G,P,'pdf',..,'levels',l) adds contours/surfaces at densities 10.^(l), with
%   default l = -2:2.
%
% KNKD_DENSITY_PLOT(G,P,'cdf',..,'quantiles',q) adds contours/surfaces at densities y that match
%   quantiles q, i.e. CDF(x < y) = q(j) for each q. Defaults are [0.5 0.75 0.9 0.95 0.99 0.999].
%
% KNKD_DENSITY_PLOT(..,'sepmdl'{ARGS}) - pass additional arguments to DIFFUSE_FRACTION, to plot
%   a separation model as a reference line.

    opt.scale = 'pdf';
    opt.dimlabels = {'k_n','k_d'};
    opt.ax = [];
    opt.lim = [0,0.8;0,0.8]';
    opt.levels = -2:2;
    opt.quantiles = [0.5 0.75 0.9 0.95 0.99 0.999];
    opt.ktgrid = 0.2:0.2:1.6;
    opt.sepmdl = {}; % {(20:10:80)','sot2'};
    
    [opt,~,isdef] = getpairedoptions(varargin,opt,'dealrest',3);
        
    w = point_weights(g{:});
    [g{:}] = ndgrid(g{:});
    assert(isequal(size(P),size(g{1})));    
    d = numel(g);
    uniform = all(P == P(1),'all');
    
    if isdef.scale && isdef.levels && ~isdef.quantiles, opt.scale = 'cdf'; end
    
   if isdef.dimlabels && d > 2
       opt.dimlabels(3:d) = arrayfun(@char,double('z')-(0:d-3),'unif',0); 
   elseif numel(opt.dimlabels) ~= d
       error('Expecting %d dimension labels');
   end
   if isdef.lim && d > 2
       opt.lim(:,3) = g{3}([1,end]);
   elseif size(opt.lim,1) ~= 2 || ~any(size(opt.lim,2) ~= [3,d])
       error('Expecting 2x%d limits');
   end

    if ~isempty(opt.ax), try axes(opt.ax); ax = opt.ax; catch, opt.ax = []; end; end
    if isempty(opt.ax)
        GUIfigure('knkd_density','',[100,100,500,400]); clf(); ax = gca();
    end
    hold(ax,'on'); box(ax,'on'); 
    axis(ax,'equal'); grid(ax,'on');
    xlabel(ax,['$' opt.dimlabels{1} '$'],'interpreter','latex','fontsize',14);
    ylabel(ax,['$' opt.dimlabels{2} '$'],'interpreter','latex','fontsize',14); 
    
    cb = colorbar();
    switch opt.scale
        case 'pdf'
            LVL = opt.levels;
            ticklbl = arrayfun(@(x) sprintf('10^{%d}',x),-LVL,'unif',0);
            set(cb.Title,'String','$pdf(x)$','Interpreter','latex','fontsize',12)        
            
        case 'cdf'
            % if isdef.levels && d > 2, opt.levels = opt.levels(1:2:end); end
            
            LVL = prctilew(P,(1-opt.quantiles)*100,P.*w);
            LVL = unique(LVL(LVL > 0));
            LVL = -log10(LVL);
            ticklbl = arrayfun(@(x) sprintf('%g%%',x),100*opt.quantiles,'unif',0);
            set(cb.Title,'String','$cdf(x)$','Interpreter','latex','fontsize',12)

        otherwise, error('Bad scale');
    end
    LVL = unique(LVL,'sorted');
    uniform = uniform || numel(LVL) < 2;
    ncolors = max(10,2*numel(LVL)+1);
    colormap(ax,parula(ncolors));
    
    if d == 2
        if ~isempty(opt.ktgrid), plotktgrid(opt.ktgrid,ax); end
        if ~isempty(opt.sepmdl), plotseparationmodel(ax,opt.sepmdl{:}); end

        if ~uniform
            [c,h] = contour(ax,g{1},g{2},-log10(P),LVL);
        end

        % h = colorbar(ax);
    else
        P = permute(P,[2,1,3:d]);
        g = cellfun(@(x) permute(x,[2,1,3:d]),g,'unif',0);
        view(-40,15);

        if d > 3
            % just pick a slice
            k = arrayfun(@(d) round(size(P,d)/2),4:d,'unif',0);            
            P = P(:,:,:,k{:});
            g = cellfun(@(x) x(:,:,:,k{:}),g(1:3),'unif',0);

            zlabel(['$' opt.dimlabels{3} '$ (slice)'],'interpreter','latex','fontsize',14);
        else
            % cpdfx = permute(cpdfx,[2,1,3]);
            zlabel(['$' opt.dimlabels{3} '$'],'interpreter','latex','fontsize',14);
        end
        if ~uniform
            for iso = LVL % isosurfaces with pdf = 0.005,0.01,0.015
                isosurface(g{:},-log10(P),iso); alpha 0.3;
            end
        end
        
        % caxis([LVL(1)-0.5,LVL(end)+0.5])
        set(gca,'OuterPosition',[0.02 0 0.8 1.0]);
        % h = colorbar();
        set(cb,'position',[0.82,0.25,0.04,0.5])
        zlim(opt.lim(:,3)');
    end
    xlim(opt.lim(:,1)');
    ylim(opt.lim(:,2)');

    if ~uniform
        caxis(ax,[LVL(1),LVL(end)]+[-0.5,0.5]*(LVL(end)-LVL(1))/ncolors);
        cb.Ticks = LVL;
        cb.TickLabels = ticklbl;
    else
        caxis(ax,[0,1]);
        cb.Ticks = [0,1];
        % cb.TickLabels = ticklbl;
        
        c = [];
        h = [];
    end

    drawnow();
end

function plotktgrid(kt,ax)

    validateattributes(kt,{'numeric'},{'positive','increasing','size',[1,NaN]'},'','ktgrid');
    x = zeros(3,numel(kt));
    x(1,:) = kt;
    x(3,:) = NaN;
    y = [flipud(x(1:2,:));x(3,:)];
    
    plot(ax,x(:),y(:),'color',[1,1,1]*0.9);
end

function plotseparationmodel(ax,varargin)
    % if nargin < 3, sunel = []; mdl = 'reindl1'; else, mdl = 'sot2'; end

    kt = 0:0.01:1.6;
    % [kd,kn] = diffuse_fraction(kt,90-sunel,mdl);
    % mdl = cellstr(mdl);
    
    % [kd,kn] = cellfun(@(m) diffuse_fraction(kt,m),mdl,'unif',0);
    % kd = kt'.*cat(1,kd{:})'; kd(end+1,:) = NaN;
    % kn = cat(1,kn{:})'; kn(end+1,:) = NaN;
    [kd,kn] = diffuse_fraction(kt,varargin{:});
    kd = kd.*kt;
    if size(kd,1) > 1
        kd(:,end+1) = NaN; kd = kd'; 
        kn(:,end+1) = NaN; kn = kn';
        
    end
    
    plot(ax,kn(:),kd(:),'color',[1 0.7 1]);
end