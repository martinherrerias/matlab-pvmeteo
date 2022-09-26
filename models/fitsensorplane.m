function [azm,tilt,calib,albedo,info] = fitsensorplane(GTI,DHI,BNI,ENI,sunel,sunaz,albedo,varargin)
% [AZM,TILT,CAL] = FITSENSORPLANE(T,LOC,GTI,GHI,DHI,BNI,SUNEL,SUNAZ,ALBEDO) - take a series of
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
    
    [arg,varargin] = getflagoptions(varargin,{'-plot'});
    
    n = size(GTI,1);
    validateattributes(GTI,{'numeric'},{'real','2d'});
    cellfun(@(x,tag) validateattributes(x,{'numeric'},{'real','vector','size',[n,1]},'',tag),...
        {DHI,BNI,ENI,sunaz,sunel},{'DHI','BNI','ENI','sunaz','sunel'});
    if knownalbedo
        validateattributes(albedo,{'numeric'},{'real','positive','vector'});
        compatiblesize(DHI,albedo);
    end
    
    filter = DHI > 10 & BNI >= -2 & sunel > 4;
    if knownalbedo && ~isscalar(albedo), filter = filter & albedo >= 0; end
        
    GTI = GTI(filter,:);
    DHI = DHI(filter);
    BNI = BNI(filter);
    ENI = ENI(filter);
    sunaz = sunaz(filter);
    sunel = sunel(filter);
    if knownalbedo && ~isscalar(albedo), albedo = albedo(filter); end
    
    Ns = size(GTI,2);

    % [F1,F2,sF1,sF2,rF1F2,Sg] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel,varargin{:});
    [F1,F2] = pvlmod_perezcoeffs(BNI,DHI,ENI,sunel,varargin{:});

    % Allow dot-product-like behaviour for arrays of any (compatible) size, without creating
    % or processing unnecessary data copies:
    sph2cartC = @(az,el) {cosd(az).*cosd(el),sind(az).*cosd(el),sind(el)};
    dotC3 = @(A,B) A{1}.*B{1} + A{2}.*B{2} + A{3}.*B{3};
    
    % S contains unit vectors pointing towards the sun
    S = sph2cartC(sunaz,sunel);
    
    % if k is a unit vector pointing towards the zenith, k'u = U{3}, k's = S{3}
    GHI = max(0,DHI + BNI.*S{3});

    ISOh = (1-F1).*DHI;                 % Isotropic component, horizontal, ISOt = ISO·(1 + k'u)/2
    CSn = F1.*DHI./max(0.087,S{3});     % Circumsolar component, normal, CSt = CSn·max(0,u's)
    HBp = F2.*DHI;                      % Horizon-brightening, vertical, HBt = HBp·sqrt(1-(k'u)²)
    ALBc = GHI;                         % Albedo, vertical, ALBt = rho·ALB·(1 - k'u)/2
    % ALBc = GHI.*(1+sind((90-sunel)/2).^2);
        
    overcast = BNI <= 0 & F1 <= 0.05;
    if ~knownalbedo, albedo = 0.2; end
    
    % Define a space over [tilt,azimuth] for each sensor
    x0 = seedvalues(GTI,S,CSn,HBp,BNI,ISOh,ALBc*albedo);
    x0 = x0(:)';
    
    lb = -Inf(1,2*Ns);
    ub = Inf(1,2*Ns);
        
    calib = ones(1,Ns);
    opt = optimset('display','off','maxfunevals',1000,'TolX',sind(0.1),'TolFun',1e-4);

    if knownalbedo
        ALBc = ALBc.*albedo;
        costfcn = @(x) GTI_error(x(1:Ns),x(Ns+1:2*Ns),1);
    else   
        x0 = [x0,0.2];
        lb = [lb,0];
        ub = [ub,1];
        costfcn = @(x) GTI_error(x(1:Ns),x(Ns+1:2*Ns),x(end));
    end
    
    XGROUPS = 50; % rough DF of Perez model... doesn't matter too much
    n = size(GTI,1);
    gidx = discretize(rand(n,1),XGROUPS);
    
    [x,info.resnorm,info.resid,info.exitflag,info.output,info.lambda,info.J] = lsqnonlin(@(x) splitapply(@std,costfcn(x),gidx),x0,lb,ub,opt);
    info.confidence_intervals = nlparci(x,info.resid,'jacobian',info.J);
    
    tilt = x(1:Ns);
    azm = x(Ns+1:2*Ns);

    badconv = tilt < 0;
    tilt(badconv) = -tilt(badconv);
    azm(badconv) = azm(badconv)+180;
    
    if ~knownalbedo, albedo = x(end); end
    
    % Refine 2nd order finite differences
    DELTA_MIN = [2,45];
    delta = min(max(diff(info.confidence_intervals(1:2*Ns,:),1,2)/2,DELTA_MIN(1)),DELTA_MIN(2));
    delta = reshape(delta,[],2)';
    
    E = zeros(3,3,Ns);
    for j = 1:3
        for k = 1:3
            z = tilt + (j-2)*delta(1,:);
            y = azm + (k-2)*delta(2,:);
            e = GTI_error(z,y,albedo);
            E(j,k,:) = std(e,1);
        end
    end
    delta(1,:) = delta(1,:)./sqrt(squeeze(diff(E(:,2,:),2,1)))';
    delta(2,:) = delta(2,:)./sqrt(squeeze(diff(E(2,:,:),2,2)))';
    delta(2,:) = min(delta(2,:),180);
    
    info.confidence_intervals = x(:) + [-1,1].*delta(:);
 
    if arg.plot
        [xx,yy] = ndgrid(linspace(-1,1,11),linspace(-1,1,11));
        E = zeros(numel(xx),Ns);
        for j = 1:numel(xx)
            x = tilt + xx(j)*delta(1,:);
            y = azm + yy(j)*delta(2,:);
            e = GTI_error(x,y,albedo);
            E(j,:) = std(e,1);
        end
        E = reshape(E,[size(xx),Ns]);

        GUIfigure('fitsensorplane');
        ax = arrayfun(@(j) subplot(1,Ns,j),1:Ns);
        for k = 1:Ns
            imagesc(ax(k),tilt(k) + xx(:,1)*delta(1,k),azm(k) + yy(1,:)*delta(2,k),E(:,:,k))
            
            xlabel(ax(k),'tilt (deg)')
            ylabel(ax(k),'azimuth (deg)')
            title(ax(k),sprintf('GTI_%d',k))
            h = colorbar(ax(k));
            h.Label.String = 'STD [W/m^2]';
        end
    end

    function E = GTI_error(surftilt,surfaz,albedo)
        U = sph2cartC(surfaz,90-surftilt);
        
        a = max(0,dotC3(S,U));       % a = max(0,u's) 
        ISO = ISOh.*(1 + U{3})/2;    %	Isotropic component: (1-F1)·DHI·(1 + k'u)/2	
        CS = CSn.*a;                 % Circumsolar component: F1·DHI·a/b
        HB = HBp.*hypot(U{1},U{2});  % Hor. brightening component: F2·DHI·sqrt(1-(k'u)²) 
        BTI = BNI.*a;                % Beam component

        ALB = albedo*ALBc.*(1 - U{3})/2;
        GTIm = ISO + CS + HB + ALB + BTI;
        
        f = overcast & isfinite(GTIm) & isfinite(GTI);
        if ~any(f), calib = 1;
        else
            calib = dot(GTIm,f)./dot(GTI,f);
        end
        
        E = double(ISO + CS + HB + ALB + BTI - GTI.*calib);
        E(~isfinite(E)) = 0;

%         A = DHI.*(a./b - (1 + U{3})/2);
%         B = DHI.*hypot(U{1},U{2});
% 
%         % Does not include covariance (several sensors at the same time)!
%         Se = sqrt(Sg.^2+(sF1.*A).^2+(sF2.*B).^2 + 2*rF1F2.*A.*B.*sF1.*sF2);
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
