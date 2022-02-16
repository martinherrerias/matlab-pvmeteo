function [MD,B,U] = bestimate(MD,varargin)
% MD = BESTIMATE(MD) - PROVISIONAL returns a constrained (closure-compliant) inverse-variance-
%   weighted estimate of kt, kd, kn (kt = kd + kn) given any number of GHI, DHI, and BNI sensors.
%
% FUTURE: Combine meteorological measurements, sensor-uncertainty information, separation- and
%   transposition-models, to deliver an Optimal State Estimate of the 'true' underlying irradiance 
%   distribution and MAP corrected sensor readings.
%
%   + Consider tilted sensors.
%   + Account for correlation between samples (sensor calibration & offset), and between sensors
%   (e.g. cosine-response, spectral-response, etc.), 
%   + Consider expected kt-kd distribution, e.g. from probabilistic separation model
%   + Account for uncertainty in closure equation for non-finite (downsampled) intervals!

    opt.discardflags = {'^','BSRN_rare_hi','BSRN_rare_lo'};
    opt = getpairedoptions(varargin,opt,'restchk');
    
    opt.discardflags = parselist(opt.discardflags,MD.flags.flags);

    MD = getsunpos(MD,0,0);
    MD = refresh(MD);
    
    if isempty(MD.uncertainty), MD = getuncertainty(MD); end
        
    % if ~all(isfield(MD,{'kt','kd','kn','F1'})) || ~all(isfield(MD.uncertainty,{'kt','kd','kn'}))
        MD = firstguess(MD,opt.discardflags);
    % end
    
    if nargout > 1
        [B,U] = unpackestimate(MD);
    end
end

function [x0,u0] = IVWmean(x,u,dim)
% Inverse-variance weighted X, according to uncertainty U

    if nargin < 3, dim = 2; end
    [x,u] = compatiblesize(x,u);
    if size(x,dim) == 1, x0 = x; u0 = u; return; end
    
    w = 1./u.^2;
    w(isnan(x)) = NaN;
    u0 = 1./sum(w,dim,'omitnan');
    x0 = sum(x.*w,dim,'omitnan').*u0;
    u0 = sqrt(u0);
end

function MD = firstguess(MD,discardflags)
% Given a set of measurements GHI, DHI, BNI with uncertainties UNC.GHI, UNC.DHI, UNC.BNI; get a 
% maximum-likelihood estimate of the 'true' G0, D0, B0. The estimates must satisfy the closure 
% condition GHI = DHI + BNI·cosZ while maximizing the log-likelyhood:
%
%       L = ln{P(G0-GHI)·P(B0 - BNI)·P(D0-DHI)}

    assert(isfield(MD,'GHI') || all(isfield(MD,{'BNI','DHI'})),'GHI source required');

    U = MD.uncertainty;
    if isempty(U)
        U = struct('GHI',1,'BNI',1,'DHI',1); % use simple averages
    end
    
    cosZ = max(0,sind(MD.data.sunel));
    if isfield(MD,'CSBNI')
        CSBNI = max(MD.data.CSBNI,[],2);
    else
        [~,CSBNI] = pvlmod_clearsky_ineichen(MD.location,MD.t,2,MD.data.sunel,[],MD.data.ENI);
    end

    % Start with inverse-variance weighted estimates
    MDc = meteoQC.flagged2nan(MD,discardflags,{'GHI','BNI','DHI'});
    
    if isfield(MD,'GHI'), [S.GHI,U.GHI] = IVWmean(MDc.data.GHI,U.GHI); end
    if isfield(MD,'DHI'), [S.DHI,U.DHI] = IVWmean(MDc.data.DHI,U.DHI); end
    if isfield(MD,'BNI'), [S.BNI,U.BNI] = IVWmean(MDc.data.BNI,U.BNI); end
    clear MDc

    if all(isfield(S,{'GHI','DHI','BNI'}))

        [S,U] = closedestimate(S,U,cosZ);
        
        % complete individual missing values from remaining two
        [S,U] = BNIfromRest(S,U,cosZ,CSBNI);
        [S,U] = GHIfromRest(S,U,cosZ);
        [S,U] = DHIfromRest(S,U,cosZ);
        
    elseif all(isfield(MD,{'GHI','DHI'})), [S,U] = BNIfromRest(S,U,cosZ,CSBNI);
    elseif all(isfield(MD,{'BNI','DHI'})), [S,U] = GHIfromRest(S,U,cosZ);
    elseif all(isfield(MD,{'BNI','GHI'})), [S,U] = DHIfromRest(S,U,cosZ);
    end

    EHI = cosZ.*MD.data.ENI;    % Extraterrestrial Horizontal Irradiance
        
    kt = S.GHI./EHI;
    kt(MD.dark) = NaN;
    kt(kt > 1.6) = NaN;

    if ~all(isfield(S,{'GHI','DHI','BNI'}))
    % Get everything from GHI, using DIRINT model for diffuse fraction

        args = {kt,'sza',90-MD.data.sunel,'Patm',MD.data.Patm}; % reindl, Orgill-Hollands, DISC, etc.
        if isfield(MD,'CSGHI')
            CSGHI = max(MD.data.CSGHI,[],2);
            ktc = CSGHI./EHI; ktc(~isfinite(CSGHI) | EHI <= 0) = NaN;
            args = [args,{'ktc',ktc,'hourangle',MD.data.hourangle}]; % Engerer2 inputs
        end
        if isfield(MD,'tpw')
            args = [args,{'tpw',MD.data.tpw,'timeseries',true}]; % DIRINT inputs
        end

        warning('off','diffuse_fraction:ignored');
        [rd,kn,mdl] = diffuse_fraction(args{:},MD.options.separationmodel);

        warning('cmd:kdguess',...
            ['Diffuse fraction estimated from GHI using the ' mdl ' model, ',...
             'this will account for significant decomposition errors.']);

        kd = rd.*kt;
        S.DHI = kd.*EHI;
        S.BNI = kn.*MD.data.ENI;
    else
        kd = S.DHI./EHI;
        kn = S.BNI./MD.data.ENI;
    end
            
    kt(kt > 1.6) = NaN;
    kn(kn > 1) = NaN;
    kd(kd > 0.95) = NaN;
    kd(kd < 0.05) = NaN;
    K = [kt,kd,kn];
    K(MD.dark,:) = NaN;
    K(MD.data.sunel < 0) = NaN;
    
    tofill = ~isfinite(K) & ~MD.dark & ~MD.missing;
    [K,~,Uk] = MeteoData.fillgaps(K,MD.timestep,hours(1),'exclude',MD.dark | MD.missing);
    K = max(K,0);

    U0 = NaN(numel(EHI),3);
    if isfield(U,'GHI'), U0 = U.GHI; end
    if isfield(U,'DHI'), U0 = U.DHI; end
    if isfield(U,'BNI'), U0 = U.BNI; end
    U0 = U0./[EHI,EHI,MD.data.ENI];
    
    Uk(~tofill) = min(U0(~tofill),Uk(~tofill));
    
    [U.kt,U.kd,U.kn] = deal(Uk(:,1),Uk(:,2),Uk(:,3));
    [S.kt,S.kd,S.kn] = deal(K(:,1),K(:,2),K(:,3));

    S.F1 = pvlmod_perezcoeffs(S.BNI,S.DHI,MD.data.ENI,MD.data.sunel);
    
    [b,MD.flags] = flagbit(MD.flags,'interp');
    tofill(:,4) = any(tofill,2);
    fld = {'kt','kd','kn','F1'};
    for j = 1:4
        MD = addfield(MD,fld{j},S.(fld{j}),['simpleguess.',fld{j}],true);
        MD.flags.(fld{j}) = bitset(MD.flags.(fld{j}),b,tofill(:,j));
    end
    
    if ~isempty(MD.uncertainty)
        for f = {'kt','kd','kn'}
            MD.uncertainty.data.(f{1}) = U.(f{1});
        end
    end
end

function [S,U] = unpackestimate(MD)

    cosZ = max(0,sind(MD.data.sunel));
    EHI = cosZ.*MD.data.ENI;
    
    S.GHI = EHI.*MD.data.kt;
    S.DHI = EHI.*MD.data.kd;
    S.BNI = MD.data.ENI.*MD.data.kn;

    if nargout > 1
        U.GHI = EHI.*MD.uncertainty.data.kt;
        U.DHI = EHI.*MD.uncertainty.data.kd;
        U.BNI = MD.data.ENI.*MD.uncertainty.data.kn;
    end
end

function [S,U] = closedestimate(S,U,cosZ,filter)
% Given a set of measurements S.GHI, S.DHI, S.BNI with uncertainties U.GHI, U.DHI, U.BNI; get a 
% maximum-likelihood estimate of the 'true' values, such that they satisfy the closure equation
% GHI = DHI + BNI·cosZ while maximizing the log-likelyhood:
%
%       L = ln{P(G0-GHI)·P(B0 - BNI)·P(D0-DHI)}
%
% Uncertainties are assumed to be independent. The solution introduces corrections on each 
% component proportional to (U.x)²/C² where U.x is the uncertainty of each component, and C the
% total uncertainty of the closure equation. 
%
% TODO: consider cos(z) uncertainty (due to "effective" solar position on finite-intervals)
    
    Ug = U.GHI.^2;
    Ud = U.DHI.^2; 
    UbC2Z = (cosZ.*U.BNI).^2;
    Uc = (Ud + UbC2Z + Ug);
    lambda = (S.GHI - S.DHI - cosZ.*S.BNI)./Uc;

    complete = isfinite(lambda);
    if nargin > 3, complete = complete & filter; end
    if ~any(complete), return; end

    lambda(~isfinite(lambda)) = 0;

    S.BNI(complete) = S.BNI(complete) + lambda(complete).*UbC2Z(complete);

    % Check and enforce B.BNI > 0, GHI >= DHI
    overcast = S.BNI < 0 & complete;
    lambda(overcast) = (S.GHI(overcast) - S.DHI(overcast))./(Ud(overcast) + Ug(overcast));
    S.BNI(overcast) = 0;

    S.GHI(complete) = S.GHI(complete) - lambda(complete).*Ug(complete);
    S.DHI(complete) = S.DHI(complete) + lambda(complete).*Ud(complete);
    
    U.GHI(complete) = U.GHI(complete).*sqrt(1 - Ug(complete)./Uc(complete));
    U.DHI(complete) = U.DHI(complete).*sqrt(1 - Ud(complete)./Uc(complete));
    U.BNI(complete) = U.BNI(complete).*sqrt(1 - UbC2Z(complete)./Uc(complete));
end

function [S,U] = BNIfromRest(S,U,cosZ,CSBNI)
% Get BNI from GHI and DHI, clip to CSBNI
    
    % ovcast = S.GHI <= S.DHI;
    % [x0,u0] = IVWmean([S.GHI(ovcast),S.DHI(ovcast)],[U.GHI(ovcast),U.DHI(ovcast)]);
    % S.GHI(ovcast) = x0;
    % S.DHI(ovcast) = x0;
    % U.GHI(ovcast) = u0;
    % U.DHI(ovcast) = u0;

    if ~isfield(S,'BNI'), missing = true(size(S.GHI)); S.BNI = zeros(numel(cosZ),1);
    else, missing = ~isfinite(S.BNI);
    end
    if ~any(missing), return; end

    % Get BNI from GHI and DHI...
    BNI = (S.GHI - S.DHI)./ cosZ;
    
    % hard clip to 0 < BNI < CSBNI
    toclip = (BNI > CSBNI);
    BNI(toclip) = CSBNI(toclip);
    
    missing = missing & isfinite(BNI);
    toclip = (toclip | BNI < 0) & missing;
    BNI = max(0,BNI);

    S.BNI(missing) = BNI(missing);
    
    U.BNI(toclip) = 0;
    [S,U] = closedestimate(S,U,cosZ,toclip);
    
    U.BNI(missing) = hypot(U.GHI(missing),U.DHI(missing))./cosZ(missing);

    % CSGHI = pvlmod_clearsky_ineichen(MD.location,MD.t,1.7,MD.sunel,MD.AMa,MD.ENI); ?
end

function [S,U] = GHIfromRest(S,U,cosZ)
% Get GHI from BNI and DHI

    if ~isfield(S,'GHI'), missing = true(size(S.BHI)); S.GHI = zeros(0,1);
    else, missing = ~isfinite(S.GHI);
    end
    if ~any(missing), return; end

    GHI = S.DHI + S.BNI.*cosZ;
    
    missing = missing & isfinite(GHI);
    GHI = max(0,GHI);

    S.GHI(missing) = GHI(missing);
    U.GHI(missing) = hypot(U.DHI(missing),U.BNI(missing).*cosZ(missing));
end

function [S,U] = DHIfromRest(S,U,cosZ)
% Get DHI from GHI and BNI, clip to DHI > 0.1·GHI

    MIN_KD = 0.1;

    if ~isfield(S,'DHI'), missing = true(size(S.GHI)); S.DHI = zeros(0,1);
    else, missing = ~isfinite(S.DHI);
    end
    if ~any(missing), return; end

    % Get BNI from GHI and DHI...
    DHI = S.GHI - S.BNI.*cosZ;
    
    toclip = DHI < S.GHI*MIN_KD;
    DHI(toclip) = S.GHI(toclip)*MIN_KD;
    
    missing = missing & isfinite(DHI);
    toclip = toclip & missing;

    S.DHI(missing) = DHI(missing);
    
    U.DHI(toclip) = 0;
    [S,U] = closedestimate(S,U,cosZ,toclip);
    
    U.DHI(missing) = hypot(U.GHI(missing),U.BNI(missing).*cosZ(missing));
end