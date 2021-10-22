function [xkd,xkn,KD,KN] = reverse_perez(surftilt,surfaz,GTI,ENI,sunel,sunaz,albedo,IAM,model,method)

% [1] D. Yang et al., “Bidirectional irradiance transposition based on the Perez model,” 
%     Solar Energy, vol. 110, pp. 768–780, Dec. 2014, doi: 10.1016/j.solener.2014.10.006.

    narginchk(6,10);
    if nargin < 7 || isempty(albedo)
        warning('Assuming 0.2 albedo');
        albedo = 0.2; 
    end
    if nargin < 8, IAM = []; end
    if nargin < 9, model = ''; end
    if nargin < 10, method = 'linear'; end
    
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
    if iscell(IAM)
        IAM = cellfun(@checkIAM,IAM,'unif',0);
        validateattributes(IAM,{'cell'},{'vector','numel',M},'','IAM');
    else
        IAM = checkIAM(IAM);
    end
    
    KDMIN = 0.05;
    NQ = 8;

    [f1c,f2c] = PerezCoefficients(model); % [8x3]

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
        e = log(e0);
        e = interpn(e,NQ)';
        e = exp(e);
        
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
    GF = {cs-iso,hb,iso + gnd,cs+gnd};
    
    A = AMr.*cz.^2.*(F1c(2,:,:).*GF{1} + F2c(2,:,:).*GF{2});
    C = cz.*((F1c(1,:,:) + z.*F1c(3,:,:)).*GF{1} + (F2c(1,:,:) + z.*F2c(3,:,:)).*GF{2} + GF{3});
    B = GF{4}.*cz;
    Gc = GTI./ENI;
    
    BC = C + B.*lambda;
    
    overcast = (A == 0);
    
    % kd = roots(A,B,C,Gc);
    
    s = BC.^2+4.*A.*Gc;
    s(s < 0) = NaN;
    s = sqrt(s);
        
    KD = cell(1,2); KN = cell(1,2);
    
    KD{1} = 0.5*(-BC + s)./A;
    KD{2} = 0.5*(-BC - s)./A;
    KD{1}(overcast) = Gc(overcast)./C(overcast);
    
    xkd = cell(1,2);
    xkn = cell(1,2);
    for j = 1:2
        KD{j}(KD{j} < KDMIN | KD{j} > 0.95) = NaN;
        KN{j} = KD{j}.*lambda;
        Kt = KD{j}+KN{j};
        [rd,Kn] = diffuse_fraction(Kt,90-sunel,'sot2');
        [xkd{j},xkn{j}] = segmentxing(KD{j},KN{j},rd.*Kt,Kn);
    end
    
    xkd = cat(3,xkd{:});
    xkn = cat(3,xkn{:});
    KD = cat(3,KD{:});
    KN = cat(3,KN{:});
    
    xkn = permute(xkn,[3,2,1]);
    xkd = permute(xkd,[3,2,1]);
    KD = permute(KD,[3,2,1]);
    KN = permute(KN,[3,2,1]);
    
    % err = forwardperez(xkd,xkn)
    
    % kd0 = mean(xkd,3,'omitnan');
    % kn0 = mean(xkn,3,'omitnan');
    % missing = ~isfinite(kd0) || ~isfinite(kn0);
    
    % if M > 1 || size(xkd) > 1 || any(missing,'all')
    %     mean(Kp./G{3},3);
    %     % lsqnonlin(
    %     (A.*kd + B).*kd - Kp
    % 
    % 
    % end
    

%     for t = 1:N
%         GUIfigure('Reverse-Perez'); clf(); hold on;
%         kt = 0:0.01:1.6;
%         [kd,kn] = diffuse_fraction(kt,90-sunel(t),'sot2');
%         kd = kd.*kt;
%         plot(kn,kd,'k:');
% 
%         plot(KN(:,:,t),KD(:,:,t));
%         plot(xkn(:,:,t),xkd(:,:,t),'o');
%     end
    
%     function err = forwardperez(kd,kn)
%         
%         kn_kd = kn./kd;
%         e = 1 + (kn./kd.*cz)./(1 + kappa.*z.^3);
%         w = {1,kd.*AMr.*cz,z};
% %         [w{:}] = compatiblesize(w{:});
% %         w = cellfun(@(x) x(:)',w,'unif',0);
% 
%         if strcmpi(method,'bin')
%             ebin = discretize(e,e0);
%             ebin(isnan(ebin)) = 1;
%             F1c = f1c(ebin,:);
%             F2c = f1c(ebin,:);
%         else
%             W = interpmatrix(e,ec,'extrap','nearest','method1d',method);
%             % W = W./sum(W,2);
%             F1c = W*f1c;
%             F2c = W*f1c;
%         end
%         F1c = reshape(F1c,[],size(e,1),size(e,2),3);
%         F2c = reshape(F2c,[],size(e,1),size(e,2),3);
%         
%         F1 = max(0,F1c(:,:,1).*w{1} + F1c(:,:,2).*w{2} + F1c(:,:,3).*w{3});
%         F2 = F2c(:,:,1).*w{1} + F2c(:,:,2).*w{2} + F2c(:,:,3).*w{3};
% 
%         err = G{1}.*F1.*kd + G{1}.*F2.*kd + G{3}.*kd + G{4}.*kn - Kp;
%     end
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

% function kd = roots(A,B,C,D)
% 
% Sbb = sum(B.^2,2);
% Sab = sum(A.*B,2);
% Sbc = sum(B.*C,2);
% Sbd = sum(B.*D,2);
% 
% a = 2*(sum(A.^2,2).*Sbb - Sab.^2);
% b = 3*(sum(A.*C,2).*Sbb - Sab.*Sbc);
% c = (sum(C.^2,2)-2*sum(A.*D,2)).*Sbb - Sbc.^2 + 2*Sab.*Sbd;
% d = Sbc.*Sbd - sum(B.*D,2).*Sbb;
% 
% kd = realcubicroots(a,b,c,d);
% 
% end

% function x = realcubicroots(a,b,c,d)
% 
%     d0 = b.^2 - 3*a.*c;
%     d1 = 2*b.^3 - 9.*a.*b.*c + 27.*d.*a.^2;
%     s = sqrt(d1.^2 - 4.*d0.^3);
%     e = ((-1 + sqrt(-3))/2).^(0:2);
%     C = e.*((d1 + sign(d1).*s)/2).^(1/3);
%     
%     realroots = abs(imag(C)) < eps(1);
%     C = real(C);
%     C(~realroots) = NaN;
%     
%     x = -(b+C+d0./C)./(3.*a);
% end


