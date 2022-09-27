function [azm,tilt,calib,albedo,info] = fitsensorplane(GTI,DHI,BNI,ENI,sunel,sunaz,albedo,varargin)
% [AZM,TILT,CAL] = FITSENSORPLANE(GTI,GHI,DHI,BNI,SUNEL,SUNAZ,ALBEDO,[TILT0,AZ0]) - take a series of
%   measurements GTI from one or several tilted radiation sensors, and perform best-fit estimates
%   of their orientation. That is, find the orientation of each sensor that minimizes the 
%   differences:
%
%           err = pvlmod_perez(tilt,azm,DHI,BNI,ENI,SUNEL,SUNAZ,ALBEDO) - cal.*GTI
%
%   Note that a calibration bias CAL accounts for possible sensitivity offsets among sensors.
%
% TODO: find also sensor fIAM, shading factors, optimize Loc.albedo, detect GHI sensor properties
% as well.
    
    if nargin < 7, albedo = []; end
    knownalbedo = ~isempty(albedo) && isfinite(albedo);
    
    % [arg,varargin] = getflagoptions(varargin,{'-plot'});
    arg.seed_tilt = [];
    arg.seed_az = [];
    arg.metric = 'rms';
    arg = parseoptions(varargin,{'-plot'},arg,'dealrest',2);
    
    arg.metric = parselist(arg.metric,{'mad','std','rms'});
    
    [n,Ns] = size(GTI);
    validateattributes(GTI,{'numeric'},{'real','2d'});
    cellfun(@(x,tag) validateattributes(x,{'numeric'},{'real','vector','size',[n,1]},'',tag),...
        {DHI,BNI,ENI,sunaz,sunel},{'DHI','BNI','ENI','sunaz','sunel'});
    if knownalbedo
        validateattributes(albedo,{'numeric'},{'real','positive','vector'});
        compatiblesize(DHI,albedo);
    end
    if ~isempty(arg.seed_tilt)
        validateattributes(arg.seed_tilt,{'numeric'},{'real','>=',0,'<=',90,'vector','numel',Ns});
    end
    if ~isempty(arg.seed_az)
        validateattributes(arg.seed_az,{'numeric'},{'real','>=',-180,'<=',180,'vector','numel',Ns});
    end
    
    filter = DHI > 10 & BNI > 0 & sunel > 4;
    if knownalbedo && ~isscalar(albedo), filter = filter & albedo >= 0; end
        
    GTI = GTI(filter,:);
    DHI = DHI(filter);
    BNI = BNI(filter);
    ENI = ENI(filter);
    sunaz = sunaz(filter);
    sunel = sunel(filter);
    if knownalbedo && ~isscalar(albedo), albedo = albedo(filter); end
    
    % [F1,F2,sF1,sF2,rF1F2,Sg] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel,varargin{:});
    [F1,F2] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel);

    % Allow dot-product-like behaviour for arrays of any (compatible) size, without creating
    % or processing unnecessary data copies:
    sph2cartC = @(az,el) {cosd(az).*cosd(el),sind(az).*cosd(el),sind(el)};
    dotC3 = @(A,B) A{1}.*B{1} + A{2}.*B{2} + A{3}.*B{3};
    
    % S contains unit vectors pointing towards the sun
    S = sph2cartC(sunaz,sunel);
    
    % if k is a unit vector pointing towards the zenith, k'u = U{3}, k's = S{3}
    GHI = max(0,DHI + BNI.*S{3});

    b = max(0.087,S{3});
    ISOh = (1-F1).*DHI;                 % Isotropic component, horizontal, ISOt = ISO·(1 + k'u)/2
    CSn = F1.*DHI./b;     % Circumsolar component, normal, CSt = CSn·max(0,u's)
    HBp = F2.*DHI;                      % Horizon-brightening, vertical, HBt = HBp·sqrt(1-(k'u)²)
    ALBc = GHI;                         % Albedo, vertical, ALBt = rho·ALB·(1 - k'u)/2
    % ALBc = GHI.*(1+sind((90-sunel)/2).^2);
        
    overcast = BNI <= 0 & F1 <= 0.05;
    if ~knownalbedo, albedo = 0.2; end
    
    % Define a space over [tilt,azimuth] for each sensor
    x0 = seedvalues(GTI,S,CSn,HBp,BNI,ISOh,ALBc*albedo);
    
    knownseeds = ~isempty(arg.seed_az) && ~isempty(arg.seed_tilt);
    if knownseeds
        c = dot(sph2cartV(x0(:,2),90-x0(:,1)),sph2cartV(arg.seed_az,arg.seed_tilt),2);
        if any(abs(c) < 0.7)
            warning('You might want to check your sensor/data azimuth convention')
        else
           x0 = [arg.seed_tilt(:),arg.seed_az(:)];
        end
    end

    x0 = x0(:)';
    
    lb = -Inf(1,2*Ns);
    ub = Inf(1,2*Ns);
        
    calib = ones(1,Ns);
    w = [];

    if knownalbedo
        ALBc = ALBc.*albedo;
        errorfcn = @(x) GTI_error(x(1:Ns),x(Ns+1:2*Ns),1);
    else   
        x0 = [x0,0.2];
        lb = [lb,0];
        ub = [ub,1];
        errorfcn = @(x) GTI_error(x(1:Ns),x(Ns+1:2*Ns),x(end));
    end
    
    switch arg.metric
    case 'mad', costfcn = @mad;
    case 'std', costfcn = @std;
    case 'rms', costfcn = @(x) sqrt(mean(x.^2));
    end
        
    opt = optimset('display','off','maxfunevals',1000,'TolX',sind(0.1),'TolFun',1e-4,'DiffMinChange',1e-1);

    [x,info.resnorm,info.resid,info.exitflag,info.output,info.lambda,info.J] = ...
            lsqnonlin(@(x) errorfcn(x),x0,lb,ub,opt);

    tilt = x(1:Ns);
    azm = x(Ns+1:2*Ns);
    
    % TODO: consult a statistician.
    % Applying weights to errors (1 = fully-directional source on the plane of the sensor)
    % and using df = sum(w) yields relatively consistent confidence bounds for the same data set 
    % at different resolutions, but is theoretically shaky.
    % The same goes for accounting for biases due to transposition.
    
        % Estimate of Perez error, due to coefficients only (standard error set to 0)
        [~,~,~,~,~,~,Se] = pvlmod_perez(tilt,azm,DHI,BNI,ENI,sunel,sunaz,albedo,[],[],[],0);

        % NLPARCI will assume errors are independent, i.e. divide by sqrt(n), but errors due to
        % transposition are likely to be consistent biases. Multiplying by sqrt(n) should account
        % for the worst-case perfect correlation.
        Se = hypot(Se.*w./vecnorm(w)*sqrt(n),info.resid).*sqrt(n./sum(w,1));

    info.confidence_intervals = nlparci(x,Se,'jacobian',info.J);

    badconv = tilt < 0;
    tilt(badconv) = -tilt(badconv);
    azm(badconv) = azm(badconv)+180;
    
    if ~knownalbedo, albedo = x(end); end

    DELTA_MIN = 2;
    if knownseeds
        c = dot(sph2cartV(azm,90-tilt),sph2cartV(arg.seed_az,90-arg.seed_tilt),2);
        DELTA_MIN = max(DELTA_MIN,acosd(min(c))*1.5);
    end
    delta = max(diff(info.confidence_intervals,1,2),DELTA_MIN);
    delta = reshape(delta,[],2)';
 
    e0 = GTI_error(tilt,azm,albedo,0);
    % c0 = costfcn(e0);
    % df = sum(w);
    if arg.plot
        colors = get(gca,'DefaultAxesColorOrder');
        
        [xx,yy] = ndgrid(linspace(-1,1,11),linspace(-1,1,11));
        E = zeros(numel(xx),Ns);
        for j = 1:numel(xx)
            x = tilt + xx(j)*delta(1,:);
            y = azm + yy(j)*delta(2,:);
            e = GTI_error(x,y,albedo);
            
            % F = c0/costfcn(e);
            % E(j,:) = 2*fcdf(F,df,df);
            E(j,:) = costfcn(e);
        end
        E = reshape(E,[size(xx),Ns]);

        GUIfigure('fitsensorplane','fitsensorplane',sprintf('%d:%d',Ns,1+knownseeds)); clf();
        if knownseeds
            ax = arrayfun(@(j) subplot(2,Ns,j),1:2*Ns);
            histopts =  {'normalization','count','displaystyle','stairs'};
        else
            ax = arrayfun(@(j) subplot(1,Ns,j),1:Ns);
        end

        for k = 1:Ns
            imagesc(ax(k),tilt(k) + xx(:,1)*delta(1,k),azm(k) + yy(1,:)*delta(2,k),E(:,:,k))
            
            hold(ax(k),'on');
            plot(ax(k),tilt(k),azm(k),'+','color',colors(1,:),'markersize',10,'linewidth',2);
            xlabel(ax(k),'tilt (deg)')
            ylabel(ax(k),'azimuth (deg)')
            title(ax(k),sprintf('GTI_%d',k))
            h = colorbar(ax(k));
            h.Label.String = [upper(arg.metric) ' [W/m^2]'];
            
            if knownseeds
                plot(ax(k),arg.seed_tilt(k),arg.seed_az(k),'+','color',colors(2,:),'markersize',10,'linewidth',2);
                e = GTI_error(arg.seed_tilt(k),arg.seed_az(k),albedo,0);
                hold(ax(k+Ns),'on');
                w0 = prctile(w(:,k),75);
                ff = w(:,k) > w0;
                ylabel(ax(k+Ns),sprintf('samples with w > %0.2f',w0));
                H = histogram(ax(k+Ns),e0(ff),histopts{:},'DisplayName','Optimized');
                histogram(ax(k+Ns),e(ff),histopts{:},'BinEdges',H.BinEdges,'DisplayName','Nominal');
                histogram(ax(k+Ns),e(ff & sunaz > 0),histopts{:},'BinEdges',H.BinEdges,'edgecolor',colors(2,:),'LineStyle',':','DisplayName','Morning');
                histogram(ax(k+Ns),e(ff & sunaz < 0),histopts{:},'BinEdges',H.BinEdges,'edgecolor',colors(2,:),'LineStyle','--','DisplayName','Evening');
                legend(ax(k+Ns),'box','off');
                xlabel(ax(k+Ns),'GTI_{mdl} - GTI_{mes} [W/m^2]');
            end
        end   
    end

    function E = GTI_error(surftilt,surfaz,albedo,tweaked)
        
        if nargin < 4, tweaked = true; end
        
        U = sph2cartC(surfaz,90-surftilt);
        
        a = max(0,dotC3(S,U));      % a = max(0,u's) 
        F_iso = (1 + U{3})/2;
        F_hb = hypot(U{1},U{2});
        ISO = ISOh.*F_iso;          %	Isotropic component: (1-F1)·DHI·(1 + k'u)/2	
        CS = CSn.*a;                % Circumsolar component: F1·DHI·a/b
        HB = HBp.*F_hb;             % Hor. brightening component: F2·DHI·sqrt(1-(k'u)²) 
        BTI = BNI.*a;               % Beam component

        ALB = albedo*ALBc.*(1 - U{3})/2;
        GTIm = ISO + CS + HB + ALB + BTI;
        
        f = overcast & isfinite(GTIm) & isfinite(GTI);
        if ~any(f), calib = 1;
        else
            calib = dot(GTIm,f)./dot(GTI,f);
        end
        
        E = double(ISO + CS + HB + ALB + BTI - GTI.*calib);
                
        if tweaked
            w = 1-dotC3(S,U).^2;    % Give more weight to points near the sensor plane
            w = w.*BTI./GTIm;       % ... with strong directional fraction
            % w = w.^2;
            E = E.*double(w);
        end
        
        switch arg.metric
        case 'mad', if tweaked, E = sign(E).*sqrt(abs(E)); end
        case 'std', E = E - mean(E);
        case 'rms'
        end

        bad = ~isfinite(E);
        if any(bad)
            E(bad) = mean(E(~bad));
        end
    end
end

function orientation = seedvalues(GTI,S,CSn,HBp,BNI,ISOh,ALB)
% Refine seed values using an approximate, linear distribution of irradiance:
% Note Horizon-brightening is placed on a point v on the ground, right below the sensor
%   GTI ~ (BNI + CSn)s·u + ISO(1+k·u)/2 + ALB(1-k·u)/2 + HBp·w·u , v = k x u x k
%   GTI - ISO/2 + ALB/2 = [(BNI + CSn)s + 1/2(ISO - ALB)k + HBp·v]·u
%   GTI - D = [B s + D k + HBp v]· u , where B = BNI + CSn, D = (ISO - ALB)/2
%   This is a constrained least-squares problem b = Ax | x'x = 1
%   Note that for points where s·u < 0, the term (BNI + kn·DHI) must be dropped

    Bs = (CSn + BNI).*[S{:}];       % Beam + CSt, spread in x,y,z components
    D = 0.5*(ISOh - ALB);  % Isotropic + albedo

    Ns = size(GTI,2);
    orientation = zeros(Ns,2); % tilt,az
 
    u0 = spherepoints(32,'regular',1)';
    u0(:,u0(3,:) < 0) = [];
    n_starts = size(u0,2);
    
    opt_rough = optimset('TolX',sind(1),'TolFun',1,'MaxFunEvals',1000);
    opt_fine = optimset('TolX',sind(0.1),'TolFun',1e-2,'MaxFunEvals',1000);

    for j = Ns:-1:1
        c = double(GTI(:,j) - D); % GTI - D
        c = max(c,0);
        dark = c == 0;
        
        ufit = NaN(3,n_starts);
        se = Inf(1,n_starts);
        its = zeros(1,n_starts+1);
        
        for k = 1:n_starts + 1
            % u = double(sph2cartV(surfaz(j),90-surftilt(j)))'; % guess value for u
            % u0 = u; % remember original guess
            if k <= n_starts
                u = u0(:,k);
                opt = opt_rough;
                TOL = 1-cosd(3);
                MAXITER = 100;
            else
                [~,i] = min(se);
                u = ufit(:,i);
                opt = opt_fine;
                TOL = 1-cosd(0.5);
                MAXITER = 100;
            end

            for i = 1:MAXITER
                prev_u = u;
                facingsun = [S{:}]*u > 0;

                % Get modeled tilted irradiance, in x,y,z components
                v = u; v(3) = 0; if norm(v) > 1e-6, v = v/norm(v); end
                A = double(D*[0,0,1] + HBp*v' + Bs.*facingsun);
                A(dark,:) = 0;

                % if iter == 1
                % On first iteration, try to improve seed with unconstrained solution
                    uu = A\c; uu = uu/norm(uu);
                    if rms(A*uu-c) < rms(A*u-c)
                        u = uu; 
                    end
                % end

                % Get constrained optimum orientation
                [u,~,exitflag] = fminsearch(@(x) double(rms(A*x/norm(x)-c)),u,opt);

                if exitflag < 1, break; end
                u = u/norm(u);
                if 1-u'*prev_u < TOL, break; end
                if i == MAXITER, exitflag = 0; end
            end
            
            its(k) = i;
            if exitflag == 1
                ufit(:,k) = u;
                se(k) = rms(A*u-c);
            end
        end
        
        orientation(j,1) = acosd(u(3));
        orientation(j,2) = atan2d(u(2),u(1));
        
%         % info.rotation(j) = acosd(u0'*u);
%         info.exitflag(j) = exitflag;
%         info.R2(j) = 1-var(A*u - c)/var(GTI(~dark,j));
%         info.MBE(j) = mean(A*u - c);
%         info.RMS(j) = rms(A*u - c);
%         m = (c + D)./(A*u + D);
%         info.slope(j) = mean(m,'omitnan');
%         v = polygon3d.rotmat([180 + S(j).Az,S(j).Tilt],'ZX')*[0;0;1];
%         scatter(A*u + D,c + D,1);
%         pause()
    end
%     disp(info) [x,resnorm,exitflag,output,lambda,J,H] = fmincon(costfcn,x0,[],[],[],[],lb,ub) 
end
