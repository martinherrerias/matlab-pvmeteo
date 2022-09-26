function varargout = plot(MD,type,varargin)
% Diagnostic plots for MeteoData objects
%
%   PLOT(MD,'ktrd',...) - kt vs rd scatter plot (colored by density)
%   PLOT(MD,'knkd'/'ktrd',...) - kn vs kd (bzw. kt vs rd) scatter plot (colored by density)
%   PLOT(MD,'shade',...) - sunel vs sunaz plot (colored by max(kn).
%   PLOT(MD,'heatmap',...) - Heat-maps for kn,kd,Ta (optionally other variables).
%
% H = PLOT(..,H) - use/return figure handle H.

    if nargin < 2, type = 'ktrd'; end
    
    TYPES = {'ktrd','kt[ -_]?.d|.d[ -_]?kt';
             'knkd','kn[ -_]?.d|.d[ -_]?kn';
             'shading','shade|shading';
             'heatmap','heat.?map'};
             % 'series','series'};
    
    type = parselist(type,TYPES,'-regexp');
        
    if ~isempty(varargin) && isscalar(varargin{end}) && ishandle(varargin{end}) && ...
       any(cellfun(@(c) isa(varargin{end},c),{'matlab.graphics.axis.Axes','matlab.ui.Figure'}))
   
        h = varargin{end};
        varargin(end) = [];
    else
        if strcmp(type,'shading'), args = {'3:1'}; else, args = {}; end
        h = GUIfigure(type,'Meteo-Data',args{:});
    end
    if nargout > 0, varargout = {h}; end
    
    switch type
        case 'heatmap', heatmaps(h,MD,varargin{:});
        case 'shading', shadingplot(h,MD,varargin{:});
        % case 'series', plotseries(h,MD,varargin{:});
        case {'ktrd','knkd'}, kxkyscatter(h,MD,type,varargin{:});
    end
    
    if getSimOption('exportplots'), exportfigure(type,h); end
end
function shadingplot(ax,MD,varargin)
    if isa(ax,'matlab.ui.Figure')
        clf(ax); ax = axes(ax);
    end
    
    opt.var = 'kn';
    opt.res = hours(MD.timestep)*15;
    opt.fcn = @max;
    opt = getpairedoptions(varargin,opt,'dealrest');
    validateattributes(opt.res,'numeric',{'scalar','real','positive','<=',15});
    validateattributes(opt.var,'char',{'nonempty'});
    
    assert(isfield(MD,opt.var),['Variable ' opt.var ' not found']);
    
    % XBINS = -180:0.4:180;
    % YBINS = 0:0.4:90;
    n = max(12,round(90/opt.res));
    XBINS = linspace(-180,180,4*n+1);
    YBINS = linspace(0,90,n+1);
    
    XTICKS = [-180:20:-140,-130:10:-50,-40:20:40,50:10:130,140:20:180];
    YTICKS = [0:10:20,30:20:90];
    
    x = (XBINS(1:end-1) + XBINS(2:end))/2;
    y = (YBINS(1:end-1) + YBINS(2:end))/2;
    
    az = solarposition.fixazimuth(MD.data.sunaz,'N2E','Eq0',MD.location.latitude);
    ix = discretize(az(~MD.dark),XBINS);
    iy = discretize(MD.data.sunel(~MD.dark),YBINS);
    z = double(MD.data.(opt.var));
    Z = accumarray([ix,iy],z(~MD.dark),[numel(XBINS),numel(YBINS)]-1,opt.fcn,NaN);

    [xx,yy] = ndgrid(x,y); 
    gaps = isnan(Z);
    if nnz(~gaps) < 3
        warning('No valid data to plot');
        G = @(x,y) NaN(size(x));
    else
        G = scatteredInterpolant(xx(~gaps),yy(~gaps),Z(~gaps),'linear','none');
    end
    
    [~,d] = solarposition.hor2eq(MD.location.latitude,xx(gaps),yy(gaps),'Eq0');
    gaps(gaps) = abs(d) < 23.5;
    Z(gaps) = G(xx(gaps),yy(gaps));

    imagesc(ax,x,y,Z','alphadata',isfinite(Z'));
    ax.YDir = 'normal';
    xticks(ax,XTICKS)
    yticks(ax,YTICKS)
    grid(ax,'on');
    ax.YMinorGrid = 'on';
    % ax.YMinorTick = 'on';
    ax.XMinorGrid = 'on';
    % ax.XMinorTick = 'on';
    
    if MD.location.latitude >=0
        xlabel(ax,'Azimuth (0° = S, 90° = E)');
    else
        xlabel(ax,'Azimuth (0° = N, 90° = W)');
    end
    set(ax,'xdir','reverse')
    ylabel('Elevation (deg)');
    
end
function kxkyscatter(ax,MD,type,varargin)
% Display kd-kt distribution of data, mark outliers in red.

    opt.flagged = false;
    opt.sepmdl = {}; % {(20:10:80)','sot2'};
    
    opt.ktlim = [0,1.3];
    opt.rdlim = [0,1];
    opt.knlim = [0,0.8];
    opt.kdlim = [0,0.8];
    opt.kngrid = 0.1:0.1:1;
    opt.ktgrid = 0.1:0.1:1.6;
    opt.rdgrid = 0:0.1:1;
    opt.kdgrid = 0:0.1:0.8;
    
    opt.auxcolor = [1 0.7 1]; % [1 0.5 0.5];
    opt.gridcolor = [1,1,1]*0.9;
    
    opt.maxpts = 50e3;
    opt.discardflags = {'^','BSRN_rare_hi','BSRN_rare_lo','CIE_out','IQR'};
    
    opt = getpairedoptions(varargin,opt,'restchk');
    
    if isa(ax,'matlab.ui.Figure')
        clf(ax); ax = axes(ax);
    end

    hold(ax,'on'); 
    box(ax,'on'); 
    axis(ax,'square');
    axis(ax,'equal'); 
    grid(ax,'on');
    
    switch type
    case 'ktrd'
        xticks(opt.ktgrid);
        yticks(opt.rdgrid);
        axis(ax,[opt.ktlim opt.rdlim]);
        xlbl = 'k_t'; 
        ylbl = 'r_d';

        % plot(ax,[-1 2],[1 1],'k:',[0 0],[-1 2],'k:');
        plotkngrid(ax,opt.kngrid,opt.ktlim,opt.gridcolor);
        axis([opt.ktlim,opt.rdlim]);

        if ~isempty(MD)
            MD = meteoQC.flagged2nan(MD,opt.discardflags,{'kt','kd'});
            x = MD.data.kt;
            y = MD.data.kd./MD.data.kt;
        end
        
        if iscell(opt.sepmdl)
           [kt,~,rd] = separationmodel(opt.sepmdl{:});
           plot(kt,rd,'color',opt.auxcolor);
        end

    case 'knkd'
        xlbl = 'k_n';
        ylbl = 'k_d';
        xticks(opt.kngrid);
        yticks(opt.kdgrid);
        axis([opt.knlim,opt.kdlim]);
        plotktgrid(ax,opt.ktgrid,opt.gridcolor);

        if ~isempty(MD)
            % MD = meteoQC.flagged2nan(MD,opt.discardflags,{'kn','kd'});
            x = MD.data.kn;
            y = MD.data.kd;
        end
        
        if iscell(opt.sepmdl)
           [kt,kn,rd] = separationmodel(opt.sepmdl{:});
           plot(kn,rd.*kt,'color',opt.auxcolor);
        end
    end
    xlabel(ax,['$' xlbl '$'],'interpreter','latex','fontsize',14);
    ylabel(ax,['$' ylbl '$'],'interpreter','latex','fontsize',14); 

    if isempty(MD), return; end
    
    flagged = MD.flags.kt > 0 | MD.flags.kn > 0 | MD.flags.kd > 0 & ~MD.dark;
    
    % Plot up to opt.maxpts, otherwise it'll hang
    f = isfinite(x) & isfinite(y) & ~MD.dark;
    n = nnz(f);
    if n > opt.maxpts
        f(f) = rand(n,1) < opt.maxpts/n;
    end
    x = x(f);
    y = y(f);   
    densityplot(x,y,3,'limits',[0 1.6 0 1.1]);
    
     if opt.flagged
        flagged = flagged(f);
        scatter(x(flagged),y(flagged),1,'r');
     end
end

function plotktgrid(ax,ktgrid,color,varargin)

    validateattributes(ktgrid,{'numeric'},{'positive','increasing','size',[1,NaN]'},'','ktgrid');
    x = zeros(3,numel(ktgrid));
    x(1,:) = ktgrid;
    x(3,:) = NaN;
    y = [flipud(x(1:2,:));x(3,:)];
    
    plot(ax,x(:),y(:),'color',color,varargin{:});
end

function plotkngrid(ax,kngrid,ktlim,color,varargin)
% Plot constant kn lines

    x = linspace(0,1,101)';
    x(end+1) = NaN;
    
    kt = (ktlim(2)-kngrid).*x + kngrid;
    rd = (kt - kngrid)./kt;
    
    plot(ax,kt(:),rd(:),'color',color,varargin{:});
end

function [kt,kn,rd] = separationmodel(varargin)

    warning_disabler = naptime('diffuse_fraction:defaultmdl'); %#ok<NASGU>
    
    kt = 0:0.01:1.6;
    [rd,kn] = diffuse_fraction(kt,varargin{:});

    if size(rd,1) > 1
        kt = compatiblesize(kt,rd);
        kt(:,end+1) = NaN; kt = kt';
        rd(:,end+1) = NaN; rd = rd'; 
        kn(:,end+1) = NaN; kn = kn';
    end
end

function heatmaps(h,MD,vars,lbl)
% seriesheatmap of kt, kd [and Ta] on subplot(3,1)

    if nargin < 3, vars = intersect({'kn','kd','Ta'},fieldnames(MD),'stable'); end
    if nargin < 4, lbl = vars; end
    
    vars = cellstr(vars);
    n = numel(vars);
    
    MD = meteoQC.flagged2nan(MD,'all');
    X = zeros(MD.Nt,n);
    MD.t.TimeZone = MD.location.TimeZone;

    idx = cellfun(@(x) ismember(vars(:),x)',MD.source,'unif',0);
    bysource = any(cat(1,idx{:}),1);
    if any(bysource)
        [X(:,bysource),~,idx] = getbysource(MD,vars(bysource)); 
        vars(bysource) = MD.Properties.VariableNames(idx);
    end
    if ~all(bysource)
        [~,idx] = parselist(vars(~bysource),fieldnames(MD));
        assert(all(arrayfun(@(j) isvector(MD.data{:,j}),idx)),...
            'Expecting source or single-column variable names');
        X(:,~bysource) = MD.data{:,idx};
    end
    
    figure(h); clf(h);
    ax = arrayfun(@(j) subplot(n,1,j),1:n);

    dayonly = ismember(vars,MeteoData.varnames({'irradiance','indices'}));
    X(MD.dark,dayonly) = NaN;
    
    for j = 1:n
        [z,~,u] = seriesheatmap(MD.t,X(:,j),'step',MD.timestep,ax(j));
        % u = u(any(Kt > 0,2));
        [a,b] = bounds(u(any(isfinite(z),2)));
        ylim([floor(a),ceil(b)]);
        cb = colorbar(ax(j)); 
        ylabel(cb,lbl{j});
    end
end
	
% function varargout = plotseries(x,varargin)
% % PLOTSERIES(Y) - Divide a time series of values Y in four horizontal sub-plots, X = 1:numel(Y)
% % PLOTSERIES(X,Y) - Do the same for paired X,Y values 
% % PLOTSERIES(...,'breaks',N) - Divide the figure in N horizontal sub-plots
% % PLOTSERIES(...,'breaks',Xb) - Divide the x-axis at explicit break-points Xb
% % PLOTSERIES(...,'name',value,...) - pass other property-value pairs to PLOT()
% %
% % h = PLOTSERIES(...) - returns an N-vector of handles, one for each sub-plot
% % [h,Xb] = PLOTSERIES(...) - returns the set of x-breaks used to divide in sub-plots
% %
% % NOTE: all subplots are left on hold-all
% 
%     [opt,varargin] = getpairedoptions(varargin,{'breaks'},{4});
% 
%     if nargin < 2 || ~isnumeric(varargin{1})
%         y = x;
%         x = 1:size(y,1);
%     else
%         y = varargin{1};
%     end
%     
%     if isscalar(opt.breaks)
%         N = opt.breaks;
%         Xb = min(x) + (max(x)-min(x))*(0:N)/N;
%     else
%         Xb = opt.breaks;
%         N = numel(Xb)-1;
%     end
%     [~,B] = histc(x,Xb);
%     B(end) = B(end-1);
% 
%     range = [min(y),max(y)];
%     range = range + [-1,1]*diff(range)/20;
%     
%     h = zeros(N,1);
%     for j = 1:N
%         h(j) = subplot(N,1,j); hold all;
%         plot(x(B == j),y(B == j,:),varargin{2:end});
%         axis([Xb(j),Xb(j+1),0,0] + axis().*[0,0,1,1]);              % fix X, keep Y range
%         if ~any(isinf(range) | isnan(range)) && range(2) > range(1)
%             axis(max([Xb(j),Xb(j+1),range],axis()));                % increase Y, if necessary
%         end
%     end
%     
%     if nargout > 0, varargout{1} = h; end
%     if nargout > 1, varargout{2} = Xb; end
% end
