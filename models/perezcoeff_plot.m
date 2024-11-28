function perezcoeff_plot(varargin)
% PEREZCOEFF_PLOT() - UI plot of Perez model anisotropy coefficients F1 and F2
% PEREZCOEFF_PLOT('-raw') - plot of raw polynomial coefficients F11, F12, .. F23
% PEREZCOEFF_PLOT(F,V,..) - add a custom set of coefficients F (with variance V) - see PEREZ_FIT

    [opt,varargin] = getflagoptions(varargin,'-raw');
    
    if nargin > 0
        validateattributes(varargin{1},{'numeric'},{'real','finite','2d','size',[8,6]},'','F');
    end
    if nargin > 1
        validateattributes(varargin{2},{'numeric'},{'real','finite','2d','size',[48,48]},'','V');
        isposdef = @(x) issymmetric(x) && all(eig(x) > -eps(1));
        assert(isposdef(varargin{2}),'Variance matrix should be symetric positive semi-definite');
    end

    if opt.raw
        GUIfigure('Perez_Coeffs_mix','Perez (1990) coefficients'); clf(); 
        plot_raw(varargin{:});
    else
        GUIfigure('Perez_Coeffs_mix','Perez (1990) coefficients'); clf(); 
        plot_F1F2(varargin{:});
    end
end

function plot_F1F2(varargin)
    %set(gcf,'DefaultAxesColorOrder',hsv(numel(MODELS)));

    TICKS_F1 = 0:0.2:1.2;
    TICKS_F2 = -1.6:0.4:0.8;
    
    ax = arrayfun(@(j) subplot(2,2,j),1:4);
    arrayfun(@formataxes,ax(1:2));

    grid(ax(1),'on');
    ylabel(ax(1),'Circumsolar fraction F_1');
    ylim(ax(1),TICKS_F1([1,end])); yticks(ax(1),TICKS_F1);

    grid(ax(2),'on')
    ylabel(ax(2),'Horizon Brightening F_2');
    ylim(ax(2),TICKS_F2([1,end])); yticks(ax(2),TICKS_F2);

    grid(ax(3),'on'); hold(ax(3),'on');
    xlabel(ax(3),'Circumsolar fraction F_1');
    ylabel(ax(3),'Horizon Brightening F_2');
    xlim(ax(3),TICKS_F1([1,end])); xticks(ax(3),TICKS_F1);
    ylim(ax(3),TICKS_F2([1,end])); yticks(ax(3),TICKS_F2);

    set(ax(4),'visible',false);

    % perez_smooth(logspace(0,1,50),0.3,45,'median','pchip');
    % drawsliders(ax{2});
    
    % TODO: alternative representation using kt-kd instead of delta-epsilon 
    % AMr = pvl_relativeairmass(zenith);
    % kz3 = 1.041*(zenith*pi/180).^3;         % kappa·z³
    % e = (1 + (1-kd)./(kd.*cosd(zenith)) + kz3)./(1+kz3); 
    % delta = kdkt.*AMr*cosz;

    if ~isempty(varargin)
        options = [{'custom','median','mean'},modelnames()];
    else
        options = [{'median','mean'},modelnames()];
    end

    h = [];
    p = uipanel(gcf,'Position',ax(4).Position);
    plotcontrols({'s','s','m','m'},...
        {'Sky brightness (\delta)','Solar zenith','Coefficients','Interp. method'},...
        {[0,0.6],[0,90],options,{'bin','nearest','linear','pchip','makima'}},...
        {0.3,45,'mean','linear'},@update,'parent',p,'position',[0.1 0.1 0.8 0.8],'-smooth');
    
    
    function update(delta,zenith,model,method,~,~)

        e0 = [1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2]';
        e0 = sort([logspace(0,1,100)';e0 - eps(e0);e0 + eps(e0)]); % max e = 10 (kd > 0.1)

        [F1{1},F2{1}] = perez_smooth(e0,delta,zenith,'min');
        [F1{2},F2{2}] = perez_smooth(e0,delta,zenith,'max');
        % [F1{3},F2{3}] = perez_smooth(e0,delta,zenith,model,method);
%         if numel(varargin) == 1
%             [F1{4},F2{4}] = perez_smooth(e0,delta,zenith,'custom',method,varargin{:});
%         elseif numel(varargin) == 2
            [F1{3},F2{3},s11,s22,~] = perez_smooth(e0,delta,zenith,model,method,varargin{:});
            F1{4} = -2.*s11 + F1{3};
            F1{5} = 2.*s11 + F1{3};
            F2{4} = -2.*s22 + F2{3};
            F2{5} = 2.*s11 + F2{3};
%         end
        
        [F1c,F2c] = cellfun(@(mdl) PerezCoefficients(mdl),modelnames(),'unif',0);
        F1c = permute(cat(3,F1c{:}),[1,3,2]);
        F2c = permute(cat(3,F2c{:}),[1,3,2]);
        F1c = max(0,F1c(:,:,1) + delta*F1c(:,:,2) + zenith*pi/180*F1c(:,:,3));
        F2c = F2c(:,:,1) + delta*F2c(:,:,2) + zenith*pi/180*F2c(:,:,3);
        
        if isempty(h)
            h{1} = plot(ax(1),e0,[F1{:}]);
            h{2} = plot(ax(2),e0,[F2{:}]);
            h{3} = scatter(ax(3),F1c(:),F2c(:),10,repmat(1:8,1,10)','marker','+');
            colorbar(ax(3));
            h{4} = plot(ax(3),[F1{3:end}],[F2{3:end}]);
            % h{3} = plot(ax(3),[F1{4}{:}],[F2{4}{:}],'+');
        else
            [h{1}.YData] = deal(F1{:});
            [h{2}.YData] = deal(F2{:});
            [h{3}.XData,h{3}.YData] = deal(F1c(:),F2c(:));
            [h{4}.XData,h{4}.YData] = deal(F1{3:end},F2{3:end});
            % [h{3}.XData] = deal(F1{4});
            % [h{3}.YData] = deal(F1{4});
        end
        [h{1}.DisplayName] = deal('min','max',model,'lo','hi');
        [h{2}.DisplayName] = deal('min','max',model,'lo','hi');
    %         for j = 1:3
    %             h{1}(j).YData = repelem(F1(:,j)',2);
    %             h{2}(j).YData = repelem(F2(:,j)',2);
    %         end 
    end
end

function MODELS = modelnames()
% Leaving 'osage1988' out! From Perez et al. (1988):
% "With the exception of Osage, data are representativeof most solar geometries"
% "the very poor results obtained withthe Osage-based model in Albuquerque and Phoenix are simply
% due to the fact that Osage, unlike the two southwestern sites, included only very few high-
% epsilon events andmean, because of the least square fitting method used to derive the coefficients,
% the resolution achieved for such events is totally unsatisfactory and allows for important 
% distortions. 
% Likewise, the small performance deterioration caused in all SNLA sites by the Albany-based model
% may be explained, in part, by the higher latitude of this site and the corresponding lack of very
% low solar zenith angle events.

    MODELS = {'allsitescomposite1990','allsitescomposite1988','sandiacomposite1988',...
             'usacomposite1988','france1988','phoenix1988','elmonte1988','albuquerque1988',...
             'capecanaveral1988','albany1988'}; % 'osage1988'
end

function plot_raw(varargin)
% PLOT_RAW() - 2x3 subplots for all Fjk coefficients
% PLOT_RAW(F,V)- 

    MODELS = fliplr(modelnames()); % plot composites on top

    x = [1, repelem([1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2],2), 10]';

    set(gcf,'DefaultAxesColorOrder',hsv(numel(MODELS)));
    ax = reshape(arrayfun(@(j) formataxes(subplot(2,3,j)),1:6),3,2)';
    for j = 1:2
        for k = 1:3
            ylabel(ax(j,k),sprintf('F_{%d%d}',j,k));
        end
    end
    
    % Plot raw coefficients
    for j = 1:numel(MODELS)
        model = MODELS{j};
        [F1c,F2c] = PerezCoefficients(model);
        style = {'linewidth',contains(MODELS{j},'composite')*2+0.5};
        plotcoeffs(F1c,F2c,style)
    end
    
    if nargin > 0
        MODELS{end+1} = 'custom';
        
        B = varargin{1};
        style = {'linewidth',3,'linestyle',':','color','k'};
        plotcoeffs(B(:,1:3),B(:,4:6),style)
    end
    
    if nargin > 1        
        % keep the diagonal elements only
        S = reshape(diag(varargin{2}).^(1/2),6,8)';
        
        style = {'linewidth',0.5,'linestyle',':','color','k'};
        
        MODELS{end+1} = 'custom +P95';
        plotcoeffs(B(:,1:3) + 2*S(:,1:3),B(:,4:6) + 2*S(:,4:6),style)
        
        MODELS{end+1} = 'custom -P95';
        plotcoeffs(B(:,1:3) - 2*S(:,1:3),B(:,4:6) - 2*S(:,4:6),style)
    end

    legend(ax(end),MODELS);
    
    function plotcoeffs(F1c,F2c,style)
        unevenstairs = @(ax,Y,varargin) plot(ax,x,repelem(Y,2,1),varargin{:});
        for i = 1:3
            unevenstairs(ax(1,i),F1c(:,i),style{:});
            unevenstairs(ax(2,i),F2c(:,i),style{:});
        end
    end
end

function ax = formataxes(ax)
    hold(ax,'on');
    set(ax,'xscale','log');
    xlabel(ax,'Sky-clearness (\epsilon)');
    xticks(ax,1:10)
end
    
function [F1,F2,s11,s22,rho] = perez_smooth(epsilon,delta,zenith,mdl,method,varargin)
% [F1,F2] = PEREZ_SMOOTH(DELTA,ZENITH,MDL,METHOD) - return circumsolar (F1) and horizon-brightening
%   (F2) coefficients for sky clearness EPSILON, sky brightness DELTA and zenith angle ZENITH, 
%   according to the coefficient set (or operator*) MDL, and (optionally) using interpolation 
%   according to METHOD instead of binning and lookup.

    MODELS = modelnames(); 
    X0 = [1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2, Inf]; % bin edges (epsilon)

    if nargin < 4, mdl = 'allsitescomposite1990'; end
    if nargin < 5, method = 'bin'; end

    zenith = zenith*pi/180;
    [epsilon,delta,zenith] = compatiblesize(epsilon,delta,zenith);
    sz = size(epsilon);
    if numel(sz) == 2 && sz(1) == 1 && sz(2) > 1
    % make sure x(filter) is [n,1]
       epsilon = epsilon';
       delta = delta';
       zenith = zenith'; 
    end
    
    filter = epsilon >= 1 & delta >= 0 & zenith <= pi/2;
    n = nnz(filter);
    epsilon = epsilon(filter);
    a = [ones(n,1),delta(filter),zenith(filter)];
    clear delta zenith

    % if ~isempty(varargin), MODELS{end+1} = 'custom'; end

    if ~any(strcmpi(MODELS,mdl)) || numel(varargin) < 2
        F1s = zeros(n,8,numel(MODELS));
        F2s = zeros(n,8,numel(MODELS));
        for j = numel(MODELS):-1:1

            model = MODELS{j};
            [F1c,F2c] = PerezCoefficients(model); % [8x3]

            F1s(:,:,j) = max(0,a*F1c');   % F11 + F12.*delta + F13.*zenith [nx8]
            F2s(:,:,j) = a*F2c';          % F21 + F22.*delta + F23.*zenith
            
            FF(:,:,j) = [F1c,F2c]; 
        end
    end
    
    if nargout > 2 
        if numel(varargin) < 2
            FF = permute(FF,[3,2,1]);  % [6xnx8]
            V = cell(1,8);
            for j = 1:8
                V{j} = cov(FF(:,:,j));
            end
            V = blkdiag(V{:});
        else
            V = varargin{2};
        end
        % V is covariance for F = [F11,F12,F13,F21,F22,F23] | bin 1, [F11 ... F23] | b2 , ..
        % Group all terms for F1jF1k into V11, F2jF2k into V22, and F1jF2k into V12
        V11 = V((0:2)' + (1:6:48),(0:2)' + (1:6:48));
        V22 = V((3:5)' + (1:6:48),(3:5)' + (1:6:48));
        V12 = V((0:2)' + (1:6:48),(3:5)' + (1:6:48)); % each [24,24]
    end

    switch lower(mdl)
    case 'custom'
        F1c = varargin{1}(:,1:3);
        F2c = varargin{1}(:,4:6);
        F1b = max(0,a*F1c');
        F2b = a*F2c';
    case 'min'
        F1b = min(F1s,[],3);
        F2b = min(F2s,[],3);
    case 'max'
        F1b = max(F1s,[],3);
        F2b = max(F2s,[],3);
    case 'median'
        F1b = median(F1s,3);
        F2b = median(F2s,3);
    case 'mean'
        F1b = mean(F1s,3);
        F2b = mean(F2s,3);
    otherwise
        [F1c,F2c] = PerezCoefficients(mdl);
        F1b = max(0,a*F1c');
        F2b = a*F2c';
    end
    
    % F1,F2 are at this point [nx8] arrays of coefficients (for the 8 bins of EPSILON)

    switch lower(method)
    case 'bin'
        ebin = discretize(epsilon,X0);
        idx = sub2ind([n,8],(1:n)',ebin);
        F1 = F1b(idx);
        F2 = F2b(idx);
        
        if nargout > 2      
            % Variance elements must be weighted by repmat([1,d,z]·[1,d,z]',8,8) 
            % weight' 
            [aTa,r,c] = xTx_triu_elements(a);
            r = (ebin - 1)*3 + r';
            c = (ebin - 1)*3 + c';
            idx = sub2ind([24,24],r,c);
            
            s11 = sqrt(dot(aTa,V11(idx),2));
            s22 = sqrt(dot(aTa,V22(idx),2));
            rho = dot(aTa,V12(idx),2)./(s11.*s22);            
        end
    otherwise
        % in = epsilon <= Xc(end);
        % F1 = F1b(:,end);
        % F2 = F2b(:,end);
        % for j = find(in)'
        %     F1(j) = interp1(Xc,F1b(j,:),epsilon(j),method);
        %     F2(j) = interp1(Xc,F2b(j,:),epsilon(j),method);
        % end
        
        %   [f, gof] = fit( (1:8)',X0(1:8)'-1,'power1','startpoint',[1e-3,3]);
        %   Xc = round(feval(f,1.5:8.5)+1,3)
        Xc = [1.015;1.091;1.295;1.708;2.425;3.552;5.204;7.505]; % bin centers
            
        ebin = interpmatrix(epsilon,Xc,'extrap','nearest','method1d',method);
        ebin = ebin./sum(ebin,2);
        
        F1 = dot(ebin,F1b,2);
        F2 = dot(ebin,F2b,2);
        
        % F1 = zeros(n,1);
        % F2 = zeros(n,1);
        % for j = 1:n
        %     F1(j) = interp1(Xc,F1b(j,:),epsilon(j),method,'extrap');
        %     F2(j) = interp1(Xc,F2b(j,:),epsilon(j),method,'extrap');
        % end
        
        if nargout > 2
        % Estimate covariance of F1,F2 based on covariance V of F 
        % If F = [F11,F12,F13,F21,F22,F23]', then [F1,F2]' = [a',0; 0,a']*F, with a' = [1,del,z]
        % Each element of the covariance of [F1 F2]' is then a'a.*Qjk, where V = [Q11,Q12;Q21;Q22]

            
            [mult,r,c] = xTx_triu_elements(ones(1,24));
            idx = sub2ind([24,24],r,c)';
            
            ra = mod(r-1,3)+1;
            ca = mod(c-1,3)+1;
            rb = floor((r-1)/3)+1;
            cb = floor((c-1)/3)+1;
            
            [~,iu,ia] = unique([ra,ca],'rows');
            aTa = a(:,ra(iu)).*a(:,ca(iu));
            
            [~,iu,ib] = unique([rb,cb],'rows');
            bTb = ebin(:,rb(iu)).*ebin(:,cb(iu));
            
            W = aTa(:,ia).*bTb(:,ib).*mult;

            s11 = sqrt(W*V11(idx)');
            s22 = sqrt(W*V22(idx)');
            rho = W*V12(idx)'./(s11.*s22);
       
        end
        
    end
    F1 = revertfilter(F1,filter);
    F2 = revertfilter(F2,filter);
end

function [xTx,r,c] = xTx_triu_elements(x)
% For a row vector x, return a list of elements of triu(x'x), with all off-diagonal elements
% multiplied by 2, so that for a symmetric matrix B, B·(x'x) = dot(B(r,c),xTx)
% For a matrix X, do the same for each row X(:,j).

    n = size(x,2);
    r = nonzeros(triu(repmat(1:n,n,1)'));
    c = nonzeros(triu(repmat(1:n,n,1)));
    xTx = x(:,r).*x(:,c).*(1 + (r ~= c)');
end
