function RES = RBE(MD,TPM,usedvar,varargin)
% R = RBE(MD,MDL,USED,...) - Recursive Bayesian Estimation test for (currently) two irradadiance
%   components: [kn,kd] given any number of sensors.
%
% INPUT: 
%   MD - MeteoData object, must include all variables used by TPM
%   TPM - Transition Probability density model(s). Table generated by KNKD_DENSITY_FIT
%   USED - cellstr. list of sensors to use as measurements, default is to use all sensors.
%
%   ..,'supervised',TRUE/FALSE - By default, if MD contains fields {kn,kd}, these will be taken
%       as reference measurements, and compared with the posterior distribution (see OUTPUT).
%   ..,'mainmdl',n/key - Use a model in the list TPM other than the first.
%   ..,'minP',1e-12 - Skip likelihood calculation for points where the prior falls below this
%           threshold. Used also to clip ignorance scores (avoid log(0)).
%   ..,'ndays',N - Reduce test to N random days in the data set.
%
%   ..,'-plot' - plot recursion steps
%   ..,'-benchmark' - calculate (and plot) deterministic benchmark
%
% OUTPUT: structure with Ntx1 vectors of error metrics (compared to "true" MD.kn, MD.kd).
%
%   RES.stats - Nt cell array of structures, each the result of GriddedPDF.stats(P) for the 
%       posterior probability distribution P.
%   RES.g - the result of recover_grids(TPM,{'kn','kd'});
%
%   On supervised mode (slow!) additionally:
%
%     RES.box - Box Ordinate Transform
%     RES.ignorance - Ignorance Scores
%     RES.energy - Energy Scores
%     RES.rmse - RMSE of deterministic benchmark
%     RES.PIT - Ntx2 Probability Integral Transform (marginal kn, and kd | kn)

    opt.ndays = Inf;
    opt.mainmdl = 1;
    opt.minP = 1e-15;
    opt.resample = 1;
    opt.supervised = all(isfield(MD,{'kn','kd'}));
    % opt.method = 'bin';
    opt = parseoptions(varargin,{'-plot','-record','-benchmark'},opt,'restchk');
    if ~opt.supervised, opt.benchmark = false; end
    
    NONO_FLAGS = {'NA','num','abs_phys','rel_phys','model','BSRN_abs_lo','BSRN_abs_hi','interp','UNC_lo','UNC_hi','shaded'};
    
    MD = getderived(MD,[TPM.conditions{:}]);
    MD = meteoQC.flagged2nan(MD,NONO_FLAGS);
    
    % Parse (composite) prior models, get all external dependencies from MD
    [mainmdl,TPM,EXT,AUX] = modeldependencies(TPM,opt.mainmdl,MD);

    g = recover_grids(TPM,{'kn','kd'});
    [kn_grid,kd_grid] = ndgrid(g{:});
    
    % start with conditional distribution on cos(z) = 0
    % dawn = GriddedPDF(g,TPM{'cosz','interpolant'}{1}.Values(:,:,1));
    
    % start with uniform distribution
    dawn = GriddedPDF(g,1);
    if opt.resample > 0, Rs = resample(dawn,opt.resample); else, Rs = dawn; end

    KNOWN = {'GTI','DHI','BNI','GHI','USW'};
    
    if nargin < 3 || (isempty(usedvar) && ~iscell(usedvar))
        usedvar = getsourceof(MD,intersect(fieldnames(MD),KNOWN));
    end
    
    % Get measurements Y
    usedvar = cellstr(usedvar);
    [Y,~,var_ix] = getbysource(MD,usedvar);
    types = MD.data.Properties.VariableNames(var_ix);
    
    Ny = size(Y,2);
    
    typ = cellfun(@(x) strcmp(types,x),KNOWN,'unif',0);
    typ = cell2struct(typ',KNOWN);
    
    [~,idx] = parselist(usedvar,{MD.sensors.ID});
    S = MD.sensors(idx);
    
    if opt.benchmark && any(typ.GHI)
    % EXPERIMENTAL: treat GHI as GTI
       typ.GTI(typ.GHI) = true;
       typ.GHI(:) = false;
    end
    
    if opt.benchmark && ~all(typ.GTI)
    % TODO!- define benchmark for other sensor types
        warning('Benchmark is still not defined for non-GTI sensors');
        opt.benchmark = false;
    end
    
    % [surftilt,surfaz,IAM,sensor_list,sensor_labels] = getsensorinfo(MD.sensors,usedvar,typ);
    
    if any(typ.GTI)
        m = [S(typ.GTI).mount];
        surftilt = [m.tilt];
        surfaz = [m.azimuth];
    else
        surftilt = [];
        surfaz = [];
    end
    sensor_labels = strjoin(usedvar,'+');
    config_name = regexprep(matlab.lang.makeValidName(sensor_labels),'[\b_]','');

    if any(cellfun(@isempty,{S.fIAM})), S = checkIAM(S); end
       
    testname = strrep(mainmdl.Properties.UserData.name,'knkd_density_','');
    testname = [testname '_' mainmdl.Properties.RowNames{1} '_' config_name];

    MD.t.TimeZone = MD.location.TimeZone;
    [days,~,id] = unique(dateshift(MD.t,'end','day'));
    % worthit = accumarray(id,any(isfinite(Y),2) | MD.dark,[],@mean) > 0.8;
    
    rng(1234);
    picks = randperm(id(end),min(id(end),opt.ndays));
        
    RES.stats = cell(MD.Nt,1);
    if opt.supervised
        RES.PIT = NaN(MD.Nt,2);
        RES.BOT = NaN(MD.Nt,1);
        RES.ignorance = NaN(MD.Nt,1);
        RES.energy = NaN(MD.Nt,1);
        RES.rmse = NaN(MD.Nt,1);
    end
    RES.g = g;

    if Ny > 1 && opt.benchmark
        eg = [{single(1:Ny)},g];
        [eg{1:3}] = ndgrid(eg{:});
        ErrInterp = griddedInterpolant(eg{:},rand(size(eg{1}),'single'));
    end
    
    if Ny > 1
        searchopts = optimset('tolx',1e-4,'tolfun',1e-2,'display','off');
    end
    
    warning_resetter = naptime({'MATLAB:contour:ConstantData',...
                                'MATLAB:legend:IgnoringExtraEntries'}); %#ok<NASGU>

    if opt.plot
        lbl = {'Prior',mainmdl.texlbl{1};
               'Likelyhood',['$ P(y_i \mid x_i),\,y_i = \{' strrep(sensor_labels,'+',', ') '\}$'];
               'Posterior','$ P(x_i \mid Y_i)$'};    
        ax = plotstep(lbl);
        MD.interval = 'e';
    end
        
    np = numel(picks);
    wb = optwaitbar(0,'RBE','UI',false);
  
    R = copy(dawn);
    P = copy(dawn);
    for k = 1:np
        d = picks(k);
        wb.update(k/np,sprintf('Day %d/%d\n',k,np),'-addtime');
        
        if opt.plot && opt.record
            videoname = ['./fig/' testname '_' config_name '_' datestr(days(d),'yyyymmdd')];
            video = VideoWriter(videoname); %#ok<TNMLP>
            video.FrameRate = 5; 
            open(video)
        end
        
        R.P = dawn.P; % Initialize Posterior (time -1)
        
        if isfield(AUX,'Kt_movmean')
            AUX.Kt_movmean.reset('cpy',AUX.Kt_marginal.get(R.P)); 
        end
        
        for t = find(id == d)'

            if MD.sunel(t) < 0, continue; end
            
            % Prior
            P.P = getprior(R,TPM,EXT(t,:),AUX);
            if all(isnan(P.P),'all'), P.P = 1; end

            % Likelihood
            lkly = ( P.P > opt.minP ); 
            Q = single(lkly);
            [Q(lkly),E] = getlikelyhood(Y(t,:),S,typ,kd_grid(lkly),kn_grid(lkly),...
                                        MD.ENI(t),MD.sunel(t),MD.sunaz(t),MD.albedo(t),...
                                        surftilt,surfaz);
            % Q = Q./sum(Q.*w,1:2);
            if any(isnan(Q),'all'), Q(:) = 1; end
            
            % Posterior
            R.P = P.P.*Q;
            
            if isfield(AUX,'Kt_movmean')
                AUX.Kt_movmean.push(AUX.Kt_marginal.get(R.P));
            end

            valid = isfinite(Y(t,:));
            ny = nnz(valid);
                        
            if opt.benchmark && ((opt.plot && ny > 0) || (ny == 1))
            % Benchmark for a single sensor: intersection of reverse-Perez and separation model
            
                [KD,KN,kd0,kn0] = reverse_perez(surftilt(valid),surfaz(valid),Y(t,valid),MD.ENI(t),...
                                  MD.sunel(t),MD.sunaz(t),MD.albedo(t),{S(valid).fIAM},[],'bin');
                
                if ny == 1 && ~isempty(kd0)
                % NOTE: take the mean if there's multiple intersections
                    RES.rmse(t) = mean(hypot(kd0 - MD.kd(t),kn0 - MD.kn(t)));
                end
            else
                kd0 = []; kn0 = []; KD = []; KN = [];
            end
            % if isempty(kd0), kd0 = NaN; end
            % if isempty(kn0), kn0 = NaN; end
                
            if isfinite(MD.kn(t)) && isfinite(MD.kd(t)) && any(isfinite(E),'all')
                if ny > 1 && opt.benchmark
                % Benchmark for multiple sensors: least RMSE from reverse-Perez
                
                    BIG = 1e4;
                    E(~isfinite(E)) = BIG;
                    ErrInterp.Values(:) = BIG;
                    ErrInterp.Values(:,lkly) = E';
                    
                    % start with min RMSE on grid
                    [~,idx] = min(sum(E.^2,2));
                    idx = find(lkly,idx);
                    x = [kn_grid(idx(end)),kd_grid(idx(end))];
                    
                    f = @(x) double(ErrInterp(1:Ny,repmat(x(1),1,Ny),repmat(x(2),1,Ny)));
                    f0 = f(double(x));
                    [x,~,residual] = lsqnonlin(f,double(x),[0,0],[1,1],searchopts);
                    if rssq(residual) > rssq(f0)
                        x = [kn_grid(idx(end)),kd_grid(idx(end))];
                    end
                    kn0 = x(1); kd0 = x(2);
                    RES.rmse(t) = hypot(kd0 - MD.kd(t),kn0 - MD.kn(t));
                end

                if opt.resample > 0
                    % Rs.P = max(0,interpn(R.grids{:},R.P,Rs.ndgrids{:},'makima'));
                    rs = interpn(R.grids{:},R.P,Rs.ndgrids{:});
                    Rs.P = rs;
                else
                    Rs = R; 
                end
                                
                RES.stats{t} = Rs.stats();
                RES.stats{t}.norm = R.lastnorm;
                
                if opt.supervised
                    RES.PIT(t,1) = Rs.PIT(MD.kn(t),1);
                    RES.PIT(t,2) = Rs.PIT(MD.kd(t),2);
                    [RES.BOT(t),px] = Rs.BOT(MD.kn(t),MD.kd(t));
                    RES.ignorance(t) = -log2(max(opt.minP,px));
                    RES.energy(t) = Rs.energyscore(MD.kn(t),MD.kd(t));
                end
            end

            if opt.plot
                plotstep(ax,P,Q,Rs,MD.kn(t),MD.kd(t),kn0,kd0,KN,KD,MD.sunel(t),MD.t(t));
                drawnow()
                    if opt.record, writeVideo(video,getframe(gcf)); end
            end
        end
        if opt.record, close(video); end
    end
end
    
function ax = plotstep(ax,P,Q,R,x0,y0,xe,ye,Xt,Yt,sunel,t)

    persistent H
    try cellfun(@delete,H(2:end)); H = H(1); end %#ok<TRYNC>
        
    if nargin < 2
    % Setup
        
        H = {};
        GUIfigure('RBE','',[100 200 1600 650]); clf();
        lbl = ax;
               
        set(gcf,'color','w');
        colormap(flipud(colormap('parula')));
        ax = arrayfun(@(j) subplot(1,3,j),1:3);
        for j = 1:3
            knkd_density_plot({[0,1],[0,1]},ones(2,2),'cdf','ax',ax(j));
            hold(ax(j),'all');
            title(ax(j),[lbl{j,1},newline(),lbl{j,2}],'interpreter','latex','fontsize',12);
            ax(j).Position = ax(j).Position*[1 0 0 0;0 1 0 0;-0.1 0 1.2 0; 0 -0.1 0 1.2];
            cb = matlab.graphics.illustration.colorbar.findColorBars(ax(j));
            cb.Position = cb.Position*diag([1 1 0.5 1]);
        end
        H{1} = text(ax(1),0,-0.2,'DD-MMM-YYYY HH:MM:SS');
        
        return;
    end  
    
    H{1}.String = datestr(t);
    H{end+1} = cdfcontour(ax(1),P,'DisplayName','AKDE');
    % [~,H{1}] = contour(ax(1),x,y,log10(P),LVL);
    
    if ~isempty(Xt)
        kt = 0:0.01:1.6;
        [kd,kn] = diffuse_fraction(kt,90-sunel,'sot2');
        kd = kd.*kt;
        H{end+1} = plot(ax(1),kn,kd,'m--','DisplayName','SOT2 sep. model');
        % legend(ax(1),[H{1:2}],'AKDE','SOT2 sep. model','edgecolor','w');
        legend(ax(1),'edgecolor','w');
    end
    
    % [~,H{3}] = contour(ax(2),x,y,log10(Q),LVL);
    H{end+1} = cdfcontour(ax(2),GriddedPDF(P.grids,Q),'DisplayName','Perez (ensemble)');
    
    if ~isempty(Xt)
        H{end+1} = plot(ax(2),Xt,Yt,'m--','DisplayName','Inv. Perez (det.)');
        if size(Xt,2) == 1
            H{end+1} = plot(ax(2),kn,kd,'--','color',[0.9 0.9 0.9],'HandleVisibility','off');
        else
            set(H{end}(2:end),'HandleVisibility','off')
        end
        % legend(ax(2),cat(1,H{3:4}),'Perez (ensemble)','Inv. Perez (det.)','edgecolor','w');
        legend(ax(2),'edgecolor','w');
    end
    
    % [~,H{5}] = contour(ax(3),x,y,log10(R),LVL);
    H{end+1} = cdfcontour(ax(3),R,'DisplayName','RBE Posterior');
    H{end+1} = plot(ax(3),x0,y0,'ro','markersize',10,'DisplayName','True (measured)');
    
    m = R.mode;
    H{end+1} = plot(ax(3),m(1),m(2),'b+','markersize',10,'DisplayName','MAP');
    
    if ~isempty(Xt)
        H{end+1} = plot(ax(3),Xt,Yt,'--','color',[0.9 0.9 0.9],'HandleVisibility','off');
        if size(Xt,2) == 1
            H{end+1} = plot(ax(3),kn,kd,'--','color',[0.9 0.9 0.9],'HandleVisibility','off');
        end
    end
    
    if ~isempty(xe)
        H{end+1} = plot(ax(3),xe,ye,'mx','markersize',10,'DisplayName','Deterministic');
        % legend(ax(3),cat(1,H{5:8}),'RBE Posterior','True (measured)','MAP','Deterministic','edgecolor','w');
    else
        % legend(ax(3),cat(1,H{5:7}),'RBE Posterior','True (measured)','MAP','edgecolor','w');
    end
    legend(ax(3),'edgecolor','w');
    
    function H = cdfcontour(ax,P,varargin)
        
        PRCTILES = [50 75 95 99];
        % cb = colorbar('peer',ax);
        cbar = matlab.graphics.illustration.colorbar.findColorBars(ax);

        lvl = prctilew(P.P,(100-PRCTILES),P.wP);
        lvl = unique(-log10(lvl(lvl > 0)),'sorted');
        if numel(lvl) < 2
            caxis(ax,[0,1]);
            cbar.Ticks = [0,1];
            H = plot(0,0,'HandleVisibility','off');
            return
        end

        ncolors = max(10,2*numel(lvl)+1);
        colormap(ax,parula(ncolors));
        
        [~,H] = contour(ax(1),P.ndgrids{:},-log10(P.P),lvl,varargin{:});

        caxis(ax,[lvl(1),lvl(end)]+[-0.5,0.5]*(lvl(end)-lvl(1))/ncolors);
        cbar.Ticks = lvl;
        cbar.TickLabels = arrayfun(@(x) sprintf('%g%%',x),PRCTILES,'unif',0);
    end
end

function [main,TPM,EXT,AUX] = modeldependencies(TPM,main,MD)
% [M,B,E,A] = modeldependencies(TPM,'main',MD) - find TPM row M matching model 'main', and if it 
%   is a composite model, return the underlying models B (otherwise B = M). If any B is conditioned
%   on external variables, return these as a table E, with field names matching the keys of B.

    % PROVISIONAL: assume zero correlation among kt distributions when estimating Kt
    % the alternative (full correlation) might need using a different PDFmovmean for sKt
    CORRELATED = false;
    AVG_WINDOW = round(1/hours(MD.timestep));

    if ~isnumeric(main)
        [~,main] = is_synonym(main,TPM.Properties.RowNames);
    end
    validateattributes(main,{'numeric'},{'scalar','integer','positive','<=',numel(TPM)});
    
    if isempty(TPM.interpolant{main}) % && contains(mdl.Properties.RowNames{main},'x')
        elements = strsplit(TPM.Properties.RowNames{main},'x');
        elements = regexprep(elements,'\((.*)\)','$1');

        [~,ie] = is_synonym(elements,TPM.Properties.RowNames);
        assert(all(ie > 0),'Combination KEY does not match precalculated elements');
    else
        ie = main;
        % [~,ie] = ismember(main,mdl.Properties.RowNames);
    end
    main = TPM(main,:);
    
    % Experimental: return objects to artificially extend the state space to kt, Kt
    AUX = struct();
    if any(ismember({'lastkt','Kt'},[TPM.conditions{ie}]))
        AUX.kt_marginal = Marginal(TPM{'lastkt','interpolant'}{1}.GridVectors,@plus);
    end
    if any(ismember({'Kt','sKt'},[TPM.conditions{ie}]))
        AUX.Kt_marginal = Marginal(TPM{'Kt','interpolant'}{1}.GridVectors,@plus);
        AUX.Kt_movmean = PDFmovmean(AUX.Kt_marginal.grids{3},AVG_WINDOW,CORRELATED);
    end
    if ismember('sKt',[TPM.conditions{ie}])
        AUX.sKt = GriddedPDF(AUX.Kt_marginal.grids{3},0);
    end
    
    TPM = TPM(ie,:);

    available = fieldnames(MD);
    
    EXT = cell(1,numel(ie));
    for k = 1:size(TPM,1)
        req = setdiff(TPM.conditions{k},{'lastkd','lastkn','lastkt','Kt','sKt'});
        if isempty(req), EXT{k} = zeros(MD.Nt,0,'single'); continue; end
        
        % Solve dependence from external variables
        [ic,iv] = ismember(req,available);
        if ~all(ic)
            error('Failed to find required %s', shortliststr(TPM.conditions{k}(~ic),'variable'));
        end
        assert(all(cellfun(@(f) isvector(MD.(f)),req)),...
            'Conditioning variables in MD are expected to be vectors');

        EXT{k} = MD.data{:,iv};
    end
    EXT = table(EXT{:},'VariableNames',TPM.Properties.RowNames');
end

function Px = getprior(last_Px,mdl,EXT,AUX)
% Integrate transition density estimate MDL weighted by last estimate LAST_PX, slicing accross
% external conditioning variables EXT.

    oths = struct();
    if isfield(AUX,'Kt_movmean')
        oths.Kt = AUX.Kt_movmean.average;
    end
    if isfield(AUX,'kt_marginal')
        oths.lastkt = AUX.kt_marginal.get(last_Px.P); 
    end
    if isfield(AUX,'sKt')
        AUX.sKt.P = oths.Kt;
        oths.sKt = sqrt(AUX.sKt.cov);
    end

    Px = 1;
    for k = 1:size(mdl,1)
        if isempty(mdl.conditions{k})
            Px = Px.*mdl.interpolant{k}.Values;
            continue;
        end

        % Solve dependence from last state, distributed with P(kn,kd) = P0
        [ic,iv] = ismember({'lastkn','lastkd'},mdl.conditions{k});
        if any(ic)
            iv = iv + 2; % position of lastkn, lastkd in interpolant dimensions
            
            % make kn,kd the last dimensions...
            d = numel(mdl.conditions{k})+2;
            if ~all(ic), iv(~ic) = d+1; d = d+1; end
            dimorder = 1:d;
            dimorder(iv) = 0;
            dimorder(1:d-2) = dimorder(dimorder ~= 0);
            dimorder(d-1:d) = iv;
            G = permute(mdl.interpolant{k}.Values,dimorder);
            % ... and reduce, weighting by last_Px
            G = sum(G.*shiftdim(last_Px.P,2-d),d-1:d);
            mdl.interpolant{k}.Values = G;
            mdl.interpolant{k}.GridVectors(iv(ic)) = [];
            mdl.conditions{k}(iv(ic)-2) = [];
        end
        
        % Solve dependency with state-derived quantities, e.g. lastkt, Kt
        fld = fieldnames(oths);
        for j = 1:numel(fld)
            [ic,iv] = ismember(fld{j},[mdl.conditions{k}]);
            if ic
                iv = iv + 2;

                % make fld{j} the last dimension...
                d = numel(mdl.conditions{k})+2;
                dimorder = [1:(iv-1),(iv+1):d,iv];
                G = permute(mdl.interpolant{k}.Values,dimorder);

                % ... and reduce, weighting by the distribution last.(fld{j})
                G = sum(G.*shiftdim(oths.(fld{j}),1-d),d);
                mdl.interpolant{k}.Values = G;
                mdl.interpolant{k}.GridVectors(iv) = [];
                mdl.conditions{k}(iv-2) = [];
            end
            if isempty(mdl.conditions{k}), break; end
        end
        if isempty(mdl.conditions{k})
            Px = Px.*mdl.interpolant{k}.Values;
            continue;
        end
        
        % ... and finally with external variables, e.g. cosz, Ksat, ...
        missing = [false,false,~isfinite(EXT{:,k})];
        if any(missing)
            dimorder = 1:numel(mdl.conditions{k})+2;
            G = mdl.interpolant{k}.Values;
            G = sum(G,dimorder(missing));
            G = permute(G,[dimorder(~missing),dimorder(missing)]);
            mdl.interpolant{k}.Values = G;
            mdl.interpolant{k}.GridVectors(missing) = [];
            mdl.conditions{k}(missing(3:end)) = [];
        end
        if isempty(mdl.conditions{k})
            Px = Px.*mdl.interpolant{k}.Values;
            continue;
        end
        
        % Generate query points for whole kn-kd grid at Nt test points
        q = mdl.interpolant{k}.GridVectors;
        [q{1:2}] = ndgrid(q{1:2});
        [s(1),s(2)] = size(q{1});
        q(3:end) = arrayfun(@(x) repmat(x,s),single(EXT{:,k}),'unif',0);
        
        % ... and evaluate PDF
        Px = Px.*mdl.interpolant{k}(q{:}); % N1 x N2 x Nt
    end
end

function [Q,E] = getlikelyhood(Y,S,typ,kd,kn,ENI,sunel,sunaz,albedo,surftilt,surfaz)
    
    G.BNI = ENI.*max(kn,0);
    G.DHI = ENI.*sind(sunel).*max(kd,0);
    G.GHI = G.DHI + G.BNI*sind(sunel);

    valid = isfinite(Y); % & isfinite(U);
    E = NaN(numel(G.BNI),numel(valid));
    
    [F1,F2,sF1,sF2,rF1F2,Se] = pvlmod_perezcoeffs(G.BNI,G.DHI,ENI,sunel);
    G.CSn = F1.*kd.*ENI;

    useful = valid & typ.GTI;
    if any(useful)
        ok = valid(typ.GTI);
        
        u = uncertainty(S(useful),Y(useful),sunaz,sunel,G.BNI+G.CSn)/1.96;
            
        [a,~,ISO,CS,HB,ALB] = pvlmod_perezgeom(surftilt(ok),surfaz(ok),sunel,sunaz);

        A = (CS - ISO);
        B = HB;
        w = sind(surftilt(ok)).^2;

        ISO = (1-F1).*G.DHI.*ISO;
        CS = F1.*G.DHI.*CS;
        HB = F2.*G.DHI.*HB;
        BTI = G.BNI.*a;
        ALB = albedo.*G.GHI.*ALB;  
        
        iam = arrayfun(@(s,a) s.fIAM(acosd(a)),S(useful)',a);
        G_mdl = ISO + HB + ALB + (CS + BTI).*iam;
        
        E(:,useful) = G_mdl - Y(useful);
        
        sF1 = sF1.*G.DHI;
        sF2 = sF2.*G.DHI;
            
        if nnz(ok) == 1
            Se = sqrt(w*Se.^2+u.^2+(sF1.*A).^2+(sF2.*B).^2 + 2*rF1F2.*A.*B.*sF1.*sF2);
            Q = normpdf(E(:,useful),0,Se);
        else
            sigma = (A.*A').*shiftdim(sF1.^2,-2) + ...
                    (B.*B').*shiftdim(sF2.^2,-2) + ...
                    (A.*B'+A'.*B).*shiftdim(sF2.*sF1.*rF1F2,-2) + ...
                    diag(w).*permute(Se.^2,[3,2,1]) + ...
                    eye(nnz(ok)).*permute(u.^2,[3,2,1]);
            Q = mvnpdf(E(:,useful),0,sigma);
        end
    else
        Q = 1;
    end

    fld = fieldnames(typ);
    for j = 1:numel(fld)
        useful = valid & typ.(fld{j});
        if ~any(useful), continue; end
        
        switch fld{j}
        case 'GTI', continue; % already handled
        case {'BNI','hBNI'}, b = G.BNI+G.CSn;            
        case 'GHI', b = G.BNI+G.CSn;
        case 'DHI', b = G.CSn;
        case 'USW', b = albedo*(G.BNI+G.CSn);
        end
        [u,iam] = uncertainty(S(useful),Y(useful),sunaz,sunel,b);
        u = u/1.96;
        
        % Correct for sensor IAM
        G_mdl = G.(fld{j}).*iam;
    
        if ~isfield(typ,fld{j}) || strcmp(fld{j},'GTI'), continue; end

        useful = valid & typ.(fld{j});
        if any(useful)
            E(:,useful) = G_mdl - Y(useful);
            Q = Q.*prod(normpdf(E(:,useful)./u),2);
        end
    end
    Q(isnan(Q)) = 0;
end

function MD = getderived(MD,conditions)
% MD = GETDERIVED(MD,CONDITIONS) calculate quantities that are not explicitly in MD, but
%   are listed in CONDITIONS. 

    if ismember('kc',conditions) && ~isfield(MD,'kc')
        kc = MD.kt.*MD.ENI.*sind(MD.sunel)./MD.CSGHI;
        MD = addsource(MD,'kc',kc,strjoin(getsourceof(MD,{'kt','CSGHI'}),'/'));
        MD.flags.data.kc = bitor(MD.flags.data.kt,MD.flags.data.CSGHI);
        MD.flags = checkfield(MD.flags,MD,@(x) x >= 0,'abs_phys',{'kc'},false);
        MD.flags = checkfield(MD.flags,MD,@(x) x < 1.6,'rel_phys',{'kc'},false);
    end
    
    if ismember('cosz',conditions) && ~isfield(MD,'cosz')
        MD = addsource(MD,'cosz',sind(MD.sunel),char(getsourceof(MD,'sunel')));
        MD.flags.data.cosz = MD.flags.data.sunel;
        MD.flags = checkfield(MD.flags,MD,@(x) x >= 0,'abs_phys',{'cosz'},false);
    end
end
