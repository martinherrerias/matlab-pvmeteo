function [KD,KN,xkd,xkn] = reverse_perez(surftilt,surfaz,GTI,ENI,sunel,sunaz,varargin)
% [KD,KN] = REVERSE_PEREZ(SURFTILT,SURFAZ,GTI,ENI,sunel,sunaz,[albedo,IAM,model,method])
%   Finds approximate piece-wise-linear solution functions to the inverse-transposition problem,
%   using the Perez et al. 1990 model.
%
% [KD,KN,xkd,xkn] = REVERSE_PEREZ(...) for a single GTI sensor, (xkd,xkn) are the intersections
%   (zero, one, or more) of the solution with a separation model (Skartveith et al. 1998). For
%   several sensors, (xkd,xkn) is a least square solution for the (possibly overdetermined)
%   inverse transposition problem.
%
%   The underlying algorithm is close to that of Yang et al. [1], although the PWL
%   implementation simplifies the search for solutions.
%
%   Run REVERSE_PEREZ() and check the code on the internal TEST function for an example. 
%
% INPUT:
%
%  SURFTILT, SURFAZ, GTI - [N,M] arrays for surface tilt (degrees), azimuth (degrees), and 
%       measured tilted irradiance corresponding to M tilted sensors at N time-steps.
%       For fixed sensors, SURFTILT and SURFAZ can be [1,M] vectors. Units for GTI and azimuth
%       convention are irrelevant, as long as they are consistent with ENI and SUNAZ.
%
%  ENI, SUNEL, SUNAZ, ALBEDO - scalars or [N,1] vectors for extraterrestrial normal irradiance,
%       solar elevation (degrees), solar azimuth (degrees), and albedo.
%
%  IAM (optional) function handle or model name recognized by CHECKIAM, or M-cell-array of models
%       (one for each GTI sensor).
%
%  MODEL (optional) - coefficient set, recognized by PEREZCOEFFICIENTS.
%
%  METHOD (optional, EXPERIMENTAL) - using 'linear', 'makima', etc. as oposed to 'bin' (default), 
%       defines an interpolation method to replace the hard binning on 'epsilon' of the original 
%       Perez et al. model. This reduces artifacts (discontinuities) in the solutions, although
%       it detaches from the canonical model implementation.
%
% OUTPUT:
%   [KD,KN] - [Q,N,M] arrays of diffuse- and direct-clearness-indices. Each column pair KD(:,j,k)
%       KN(:,j,k) represents a piece-wise-linear-approximation to the solution space of the 
%       inverse transposition problem for sensor k at timestep j. That is, any point (kd,kn) 
%       interpolated from the PWL curve KD(:,j,k), KN(:,j,k)
%       
% pvl_perez(SurfTilt, SurfAz, DHI, DNI, HExtra, SunZen, SunAz, AM, varargin)
%
% [1] D. Yang et al., “Bidirectional irradiance transposition based on the Perez model,” 
%   Solar Energy, vol. 110, pp. 768–780, Dec. 2014, doi: 10.1016/j.solener.2014.10.006.
% [2] Perez, R., Ineichen, P., Seals, R., Michalsky, J., Stewart, R., 1990. "Modeling daylight 
%   availability and irradiance components from direct and global irradiance".
%   Solar Energy 44 (5), 271–289.
%
% See also: PVL_PEREZ, PVLMOD_PEREZ

    if nargin == 0, test(); return; end
            
    narginchk(6,Inf);

    opt.albedo = 0.2;
    opt.IAM = 'none';
    opt.model = [];
    opt.method = 'bin';
    
    opt.Nq = 8;
    opt.kdmin = 0.05;
    
    opt = getpairedoptions(varargin,opt,'dealrest',4);
    
    [albedo,method,NQ,KDMIN] = deal(opt.albedo,opt.method,opt.Nq,opt.kdmin);
    assert(KDMIN > 0,'kdlim(1) > 0 is required to set a max. bound on epsilon');

    validateattributes(ENI,{'numeric'},{'vector','real','size',[NaN,1]});
    validateattributes(sunel,{'numeric'},{'vector','real','size',[NaN,1]});
    validateattributes(sunaz,{'numeric'},{'vector','real','size',[NaN,1]});
    validateattributes(albedo,{'numeric'},{'vector','real','size',[NaN,1]});
    [sunel,sunaz] = compatiblesize(sunel,sunaz,ENI,albedo);
    N = numel(sunel);
    
    validateattributes(surftilt,{'numeric'},{'real','2d'});
    validateattributes(surfaz,{'numeric'},{'real','2d'});
    s = compatiblesize(surftilt,surfaz,GTI,sunel,'-size');
    validateattributes(GTI,{'numeric'},{'real','2d','size',s});
    M = s(2);
    
    warning_resetter = naptime('checkIAM:weird'); %#ok<NASGU>
    if iscell(opt.IAM)
        IAM = cellfun(@checkIAM,opt.IAM,'unif',0);
        validateattributes(IAM,{'cell'},{'vector','numel',M},'','IAM');
    elseif isempty(opt.IAM)
        IAM = checkIAM('none');
    else
        IAM = checkIAM(opt.IAM);
    end

    [f1c,f2c] = PerezCoefficients(opt.model); % [8x3]

    kappa = 1.041;
    z = (90-sunel)*pi/180;
    AMr = pvl_relativeairmass(max(0,90-sunel));

    e0 = [1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2, 1/KDMIN]';  % bin edges (epsilon)
    ec = [1.015;1.091;1.295;1.708;2.425;3.552;5.204;7.505];     % bin 'pivots' for interpolation
    
    if strcmpi(method,'bin')
        e = log([e0(1:8),e0(2:9)]);
        e = interpn(e,repmat((1:8)',1,NQ),repmat(linspace(1,2,NQ),8,1));
        e = exp(e);
        F1c = repelem(f1c,NQ,1);
        F2c = repelem(f2c,NQ,1);
    else
        e = exp(interpn(log(e0),log2(NQ)))';
        W = interpmatrix(e,ec,'method1d',method,'extrap','nearest');
        F1c = W*f1c;
        F2c = W*f2c;
    end
    
    e = reshape(e',1,1,[]);
    F1c = permute(F1c,[2,3,1]);
    F2c = permute(F2c,[2,3,1]);

    cz = max(0,sind(sunel));
    
    % kn/kd = BHI/DHI = "lambda" in Yang et al. 2014
    lambda = (e-1).*(1 + kappa.*z.^3).*cz; % bin edges (kn/kd)

    [a,b,iso,~,hb,gnd] = pvlmod_perezgeom(surftilt,surfaz,sunel,sunaz);
    
    if iscell(IAM)
        for j = 1:M
           a(:,j) = a(:,j).*IAM{j}(acosd(a(:,j))); 
        end
    else
        a = a.*IAM(acosd(a));
    end
    gnd = albedo.*gnd;
    cs = a./b;
    
    % Geometric weights for circumsolar, horizon-brightening, and isotropic components
    GF = {cs-iso, hb, iso+gnd};
    
    A = AMr.*cz.^2.*(F1c(2,:,:).*GF{1} + F2c(2,:,:).*GF{2});
    C = cz.*((F1c(1,:,:) + z.*F1c(3,:,:)).*GF{1} + (F2c(1,:,:) + z.*F2c(3,:,:)).*GF{2} + GF{3});
    B = a + gnd.*cz;
    Gc = GTI./ENI;
    
    BC = C + B.*lambda;
    overcast = (A == 0);
        
    s = BC.^2+4.*A.*Gc;
    s(s < 0) = NaN;
    s = sqrt(s);
        
    KD = cell(1,2); KN = cell(1,2);
    
    KD{1} = 0.5*(-BC + s)./A;
    KD{2} = 0.5*(-BC - s)./A;
    KD{1}(overcast) = Gc(overcast)./C(overcast);
    
    for j = 1:2
        KD{j}(KD{j} < KDMIN | KD{j} > 0.95) = NaN;
        KN{j} = KD{j}.*lambda;
    end

    if M == 1 && nargout > 2
        xkd = cell(1,2);
        xkn = cell(1,2);
    
        for j = 1:2
            Kt = KD{j}+KN{j};
            [rd,Kn] = diffuse_fraction(Kt,90-sunel,'sot2');
            [xkd{j},xkn{j}] = segmentxing(KD{j},KN{j},rd.*Kt,Kn);
        end
        
        xkd = cat(3,xkd{:});
        xkn = cat(3,xkn{:});
        xkn = permute(xkn,[3,2,1]);
        xkd = permute(xkd,[3,2,1]);
    end

    KD = cat(3,KD{:});
    KN = cat(3,KN{:});
    
    if M > 1 && nargout > 2
        G_mdl = forwardmodel(mean(KD,2,'omitnan'),mean(KN,2,'omitnan')).*ENI;

        xkd = zeros(N,1);
        xkn = zeros(N,1);
        E = mean((G_mdl - GTI).^2,2);
        [~,idx] = min(E,[],3);
        for j = 1:N
            xkd(j) = mean(KD(j,:,idx(j)),'omitnan');
            xkn(j) = mean(KN(j,:,idx(j)),'omitnan');
        end
        
        x0 = double([xkd;xkn]);
        % x = lsqnonlin(@(x) forwardmodel(x(1:N),x(N+1:end)).*ENI - GTI,x0,zeros(2*N,1),ones(2*N,1));
        x = fminsearch(@(x) mean((forwardmodel(x(1:N),x(N+1:end)).*ENI - GTI).^2,1:2),x0);
        xkd = x(1:N);
        xkn = x(N+1:end);
    end
    
    KD = permute(KD,[3,2,1]);
    KN = permute(KN,[3,2,1]);
    
    
    function kp = forwardmodel(kd,kn)
        
        e = 1 + kn./(kd.*cz)./(1 + kappa.*z.^3);
        w = {1,kd.*AMr.*cz,z};

        if strcmpi(method,'bin')
            ebin = discretize(e,[e0(1:end-1);Inf]); 
            ebin(isnan(ebin)) = 1;
            F1c = f1c(ebin,:);
            F2c = f2c(ebin,:);
        else
            W = interpmatrix(e(:),ec,'extrap','nearest','method1d',method);
            % W = W./sum(W,2);
            F1c = W*f1c;
            F2c = W*f2c;
        end
        F1c = reshape(F1c,[size(e),3]);
        F2c = reshape(F2c,[size(e),3]);
        d = ndims(F1c);
        
        [w{:}] = compatiblesize(w{:});
        W = cat(d,w{:});
        
        F1 = max(0,dot(F1c,W,d)); % F1c(:,:,1).*w{1} + F1c(:,:,2).*w{2} + F1c(:,:,3).*w{3});
        F2 = dot(F2c,W,d); % F2 = F2c(:,:,1).*w{1} + F2c(:,:,2).*w{2} + F2c(:,:,3).*w{3};
        
        kp = (((1-F1).*iso + F2.*hb + F1.*cs).*kd + (kd+kn).*gnd).*cz + a.*kn;  
    end
end

function [X,Y] = segmentxing(xa,ya,xb,yb)
% Find all intersections within the NxM pairs of piece-wise-linear curves:
%
%   PWL{xa(i,j,:),ya(i,j,:)} = PWL{xb(i,j,:),yb(i,j,:)}
%
% Return NxMxC arrays X,Y, where C is given by the max. number of intersections 
% found for any i,j.

    y34 = yb(:,:,1:end-1)-yb(:,:,2:end);
    x34 = xb(:,:,1:end-1)-xb(:,:,2:end);
    x13 = xa(:,:,1:end-1)-xb(:,:,1:end-1);
    y13 = ya(:,:,1:end-1)-yb(:,:,1:end-1);
    x1 = xa(:,:,1:end-1);
    y1 = ya(:,:,1:end-1);
    x12 = x1-xa(:,:,2:end);
    y12 = y1-ya(:,:,2:end);

    D = (x12.*y34 - y12.*x34);
    t = (x13.*y34 - y13.*x34)./D;
    u = (x13.*y12 - y13.*x12)./D;

    xing = (t >= 0) & (t <= 1) & (u >= 0) & (u <= 1);
    sz = size(xing);

    X = NaN(sz);
    Y = NaN(sz);
    X(xing) = x1(xing) - t(xing).*x12(xing);
    Y(xing) = y1(xing) - t(xing).*y12(xing);

    [xing,idx] = sort(xing,3,'descend');
    [i,j,~] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
    idx = sub2ind(sz,i,j,idx);
    X = X(idx);
    Y = Y(idx);
    
    used = any(xing,1:2);
    X = X(:,:,used);
    Y = Y(:,:,used);
end

function test()
    
    METHOD = 'bin';

    kt = rand(2,1)+0.2;
    sunel = rand(2,1)*90;
    sunaz = rand(2,1)*180 - 90;
    ENI = 1360;
    
    [rd,kn] = diffuse_fraction(kt,90-sunel,'sot2');
    DHI = rd.*kt.*ENI.*sind(sunel).*(1+rand(size(rd))*0.1);
    BNI = kn.*ENI.*(1+rand(size(kn))*0.1);

    surftilt = rand(1,3)*90;
    surfaz = rand(1,3)*180 - 90;

    % GTI = [300,600,900;800,1000,900];
    albedo = 0.2;
    IAM = checkIAM('martin-ruiz','maxloss',0.1);

    GTI = pvlmod_perez(surftilt,surfaz,DHI,BNI,ENI,sunel,sunaz,albedo);
    GTI = GTI + randn(size(GTI))*10;
    
    % rng(123);
    % reverse_perez(60,0,800,1360,45,0,0.2,IAM,[],'bin');
    [KD,KN,xkd,xkn] = reverse_perez(surftilt,surfaz,GTI,ENI,sunel,sunaz,albedo,IAM,[],METHOD);
    
    h = GUIfigure('reverse_perez_test'); clf(h);
    plot(MeteoData(),'knkd',h); hold on;
    for j = 1:size(KD,2)
        for k = 1:size(KD,3)
            plot(KN(:,j,k),KD(:,j,k)); 
        end
    end
    
    plot(xkn,xkd,'o');
end
