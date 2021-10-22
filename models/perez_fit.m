function [B,se,VB,ec] = perez_fit(GTI,surftilt,surfaz,DHI,BNI,ENI,sunel,sunaz,albedo,varargin)
% [B,se,VB] = PEREZ_FIT(GTI,surftilt,surfaz,DHI,BNI,ENI,sunel,sunaz,albedo,[IAM,method])

% TODO: adjust diffuse components & horizontal measurements for IAM
% TODO: weight observations based on (correlated?!) sensor uncertainty        

    opt = getpairedoptions(varargin,{'fIAM','method'},{@(x) x < 90,'bin'},'dealrest');
    parsestruct(opt,'method','class','char',@(x) ismember(x,{'bin','linear','pchip','makima'}));

    validateattributes(GTI,{'numeric'},{'real','2d'});
    [n,m] = size(GTI);
    cellfun(@(x,tag) validateattributes(x,{'numeric'},{'real','vector','size',[n,1]},'',tag),...
        {DHI,BNI,ENI,sunaz,sunel},{'DHI','BNI','ENI','sunaz','sunel'});
    validateattributes(albedo,{'numeric'},{'real','vector'});
    albedo = compatiblesize(albedo,DHI);
    validateattributes(surftilt,{'numeric'},{'real','vector','>=',-90,'<=',90,'size',[NaN,m]});
    validateattributes(surfaz,{'numeric'},{'real','vector','>=',-360,'<=',360,'size',[NaN,m]});
    [surftilt,surfaz] = compatiblesize(surftilt,surfaz,GTI);
    
    if nargin < 10 || isempty(opt.fIAM), fIAM = @(x) x < 90; 
    else
        warning_resetter = naptime('checkIAM:weird'); %#ok<NASGU>
        if iscell(opt.fIAM)
            fIAM = cellfun(@checkIAM,opt.fIAM,'unif',0);
            validateattributes(fIAM,{'cell'},{'vector','numel',m},'','fIAM');
        else
            fIAM = checkIAM(opt.fIAM);
        end
    end

    valid = DHI > 20 & BNI >= 0 & sunel > 4 & albedo >= 0 & ENI > 0;
    DHI = double(DHI(valid)); 
    BNI = double(BNI(valid));
    ENI = double(ENI(valid));
    GTI = double(GTI(valid,:));
    sunel = double(sunel(valid));
    sunaz = double(sunaz(valid));
    albedo = double(albedo(valid));
    if size(surftilt,1) > 1
        surftilt = double(surftilt(valid,:));
        surfaz = double(surfaz(valid,:));
    end

    AMr = pvl_relativeairmass(90-sunel);
    kappa = 1.041;
    z = max(0,90-sunel)*pi/180;
    
    e = 1 + (BNI./DHI)./(1+kappa.*z.^3);
    e(DHI == 0) = 1;
    del = DHI.*AMr./ENI;

    % if k is a unit vector pointing towards the zenith, k'u = U{3}, k's = S{3}
    GHI = max(0,DHI + BNI.*sind(sunel));

    [a,b,ISO,~,HB,ALB] = pvlmod_perezgeom(surftilt,surfaz,sunel,sunaz);

    if iscell(fIAM)
        for j = 1:m
            [IAMiso,IAMgnd,IAMhb] = diffuseIAM(fIAM{j},surftilt(:,j));
            a(:,j) = a(:,j).*fIAM{j}(acosd(a(:,j)));
            ISO(:,j) = ISO(:,j).*IAMiso;
            ALB(:,j) = ALB(:,j).*IAMgnd;
            HB(:,j) = HB(:,j).*IAMhb;
        end
    else
        [IAMiso,IAMgnd,IAMhb] = diffuseIAM(fIAM,surftilt);
        a = a.*fIAM(acosd(a));
        ISO = ISO.*IAMiso;
        ALB = ALB.*IAMgnd;
        HB = HB.*IAMhb;
    end
    ALB = albedo.*ALB.*GHI;
    ISO = ISO.*DHI;

    A{1} = DHI.*a./b - ISO;
    A{2} = A{1}.*del;
    A{3} = A{1}.*z;
    A{4} = DHI.*HB;
    A{5} = A{4}.*del;
    A{6} = A{4}.*z;
    A = cat(3,A{:});
    
        % Temps & Coulson
        % GND = 0.5*albedo.*GHI.*(1 - U{3}).*(1+sin(z/2).^2).*abs(cosd(surfaz-sunaz));
    
    % Beam tilted component
    % BTI = BNI.*a.*fIAM(acosd(a));
    BTI = BNI.*a;
    
    Y = GTI - BTI - ALB - ISO; % Y = DTI - isotropic/(1-F1)
    
    % A = A./DHI;
    % Y = Y./DHI;
    
    switch opt.method
    case 'bin'
        % Select which bin e falls into (simplified from pvl_perez)
        ebin = discretize(e,[1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2, Inf]);
        ec = accumarray(ebin,e,[],@mean);

        B = nan(6,8);
        se = nan(1,8);
        VB = cell(1,8);
        for j = 1:8
            inj = (ebin == j);
            X = reshape(A(inj,:,:),[],6);
            y = reshape(Y(inj,:),[],1);
            valid = all(isfinite(X),2) & isfinite(y);
            X = X(valid,:);
            y = y(valid);

            n = nnz(inj);
            B(:,j) = X\y;
            ey = y - X*B(:,j);
            se(j) = sqrt((ey'*ey)/(n-6));
            VB{j} = inv(X'*X).*se(j)^2;

            % [B(:,j),stats] = robustfit(X,y,[],[],0);
            % se(j) = stats.ols_s;
            % v = repmat(DHI(inj),size(Y,2),1);
            % se(j) = std(stats.resid.*v(valid),0).*sqrt(n/(n-6));
            % VB{j} = stats.covb;
        end
        VB = blkdiag(VB{:});
    otherwise

        % X0 = [1, 1.065, 1.23, 1.5, 1.95, 2.8, 4.5, 6.2, 15]; % bin edges (epsilon)
        % [f, gof] = fit( (1:8)',X0'-1,'power1','startpoint',[1e-3,3]);
        % Xc = round(feval(f,1.5:8.5)+1,3)
        Xc = [1.015;1.091;1.295;1.708;2.425;3.552;5.204;7.505]; % bin centers
            
        ebin = interpmatrix(e,Xc,'extrap','nearest','method1d',opt.method);
        ebin = ebin./sum(ebin,2);
        
        ebin = repmat(ebin,size(A,2),1);

        X = reshape(A,[],6);
        y = Y(:);
        valid = all(isfinite(X),2) & isfinite(y);
        X = X(valid,:);
        y = y(valid);
        ebin = ebin(valid,:);
        ec = sum(ebin,1);

        X = repelem(ebin,1,6).*repmat(X,1,8);

        B = nan(6,8);
         % n = numel(y);
        B(:) = X\y;
        ey = (y - X*B(:)).*ebin;
        se2 = (ey'*ey)./(ec-6);
                
        [B(:),~,~,VB] = lscov(X,y,1./(ebin*diag(se2)));
            
        % n = numel(y);
        B(:) = X\y;
        ey = (y - X*B(:)).*ebin;
        % v = repmat(DHI,size(Y,2),1);
        % ey = (y - X*B(:)).*ebin.*v(valid);
        se2 = (ey'*ey)./(ec-6);
        se = sqrt(full(diag(se2)));
        
        % VB = inv(X'*X).*se^2;
        
%         [B(:),stats] = robustfit(full(X),y,[],[],0);
%         se = smatlab copy upper triangular to lowertats.s;
%         VB = stats.covb;
    end
    