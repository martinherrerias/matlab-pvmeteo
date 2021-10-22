function MD = fitclearsky(MD,varargin)
% MD = FITCLEARSKY(MD) - Iterative detection of clear-sky samples, and of a smooth 
%   Linke-Turbidity time series that fits those detected samples. The main goal of the algorithm
%   is not to provide the closest possible clear-sky model, but to provide a model that:
%
%       a) Provides a reasonable fit to measured GHI & BNI
%       b) Can be used to estimate robust confidence intervals for BNI and DHI/GHI
%       c) Can be used to estimate effective solar positions, such that cos(Z) = CSBHI/CSBNI
%
%   While PVLMOD_CLEARSKY_INEICHEN using look-up-table turbidity-factors often fails on the first
%   points, CAMS McClear data is problematic on the last two (errors tend to be small for most of
%   the data, yet unpredictably high for particular days; and the use of no refraction correction
%   makes it problematic to use CSBHI/CSBNI for samples with low altitude).
%
%   In general terms, the algorithm is the following:
%
%     0. Start with TL = TL0 = LINKETURBIDITY(LOC,T) - i.e. Linke-Turbidity values interpolated 
%        from a look-up table. If clear-sky values (CSGHI,CSBNI) are already included in MD, 
%        calculate their equivalent Linke-Turbidity, and average with TL.
%     1. Detect clear-samples CS, using PVLMOD_DETECT_CLEAR_TIMES based on GHI and on CSGHI 
%        calculated from PVLMOD_CLEARSKY_INEICHEN(..,TL,...);
%     2. Get TL(CS), such that  CSGHI(CS) == PVLMOD_CLEARSKY_INEICHEN(..,TL(CS),...). If BNI
%        is available, calculate its equivalent Linke-Turbidity, and average with TL.
%        Remove any non-physical and outlier values of TL. Treat all TL(~CS) as missing.
%     3. Use a moving-average filter to recover a smooth TL function, i.e. fill the gaps TL(~CS)
%        while reducing its variability.
%     4. Use the new-found TL to update the estimate for CSGHI. Iterate 1-3 until convergence.
%        If the fraction of clear-sky points nnz(CS)/nnz(CSGHI > 0) is less than a threshold
%        X (see below), or if the algorithm fails to converge, TL will be set back to seed TL0.
%     5. Use the final TL to (possibly) improve clear-sky estimates by evaluating the clear-sky
%        model at 1-min resolution. See EFFECTIVESOLARPOSITION.
%
% INPUT:
%   LOC - location structure {latitude,longitude,altitude} as required by PVL_EPHEMERIS
%   T - vector of DATENUM timesteps, expected to be regular
%   MD - structure with Global Horizontal Irradiance vector GHI, and optional:
%       BNI (recommended) - Direct Normal Irradiance (used to correct TL, if available)
%       MD.Patm, MD.Ta (recommended) vectors of atmospheric pressure (Pa) and temperature (°C),
%           to account for refraction near the horizon. Defaults are 1.01e5 Pa and 10°C.
%       MD.data.CSGHI,MD.data.CSBNI (optional) used only as seed values for Linke-Turbidity.
%       MD.data.AMa, MD.data.ENI (optional) absolute air-mass and Extraterrestrial Normal Irradiance 
%   TL - (optional) vector of Linke-Turbidity values, for PVLMOD_CLEARSKY_INEICHEN
%
% OUTPUT: 
%   C - structure of results, with fields {TL,CS,CSGHI,CSBNI,CSDHI,ENI}:
%       R.CS - boolean vector of detected clear samples.
%       R.TL - vector of smoothed Linke-Turbidity values, fitted to clear-samples R.CS or reset
%           to the seed TL0, if not enough clear samples were detected.
%       R.CSGHI,R.CSBI,R.CSDHI- vectors of  clear-sky-global, -beam, and -diffuse irradiances,
%           calculated from PVLMOD_CLEARSKY_INEICHEN(LOC,T,R.TL,..).
%       R.ENI - Extraterrestrial Normal Irradiance, byproduct of EFFECTIVESOLARPOSITION.
%   SP - structure with fields {Az,El,w,dec} - vectors of effective azimuth (N2E), effective
%       apparent elevation,  hour angle, and declination (all in degrees), considering the finite
%       length of intervals [T - dT/2, T + dT/2]. See EFFECTIVESOLARPOSITION.
%
% [C,SP] = FITCLEARSKY(...,'Name',Value) - optional arguments to control the algorithm:
%
%   'nmin',N: ensure that at least N weighted-samples* are used for averaging any given timestep. 
%       If the density of clear-samples at any point does not meet this criterion, the averaging
%       window size for the whole series will be increased (to get an estimate or uncertainty).
%       (*) NOTE: Points are weighted by W = 4/max(2,AM)², so that timesteps with 1 < AM < 2 count
%       as 1 sample, but e.g. points with AM = 4 count as 1/4 of a sample.
%       Default (recommended) is N ~ (24 weighted-hours / dT).
%
%   'P',P: return TL = TL* - PRCTILE(TL* - TL0,P), where TL* is the smooth moving-average, 
%       and TL0 is the calculated clear-sample turbidity. 
%       By default P = 50, i.e. TL is adjusted so that the median of the errors TL* - TL0 = 0. 
%       Higher values (e.g. 90, 95) can be used to get minimum expected values for turbidity.
%
%   'minGHI',20 [W/m²],'minsunel',8 [°],'limTL',5.0 - Modify the filters used to select meaning-
%       full points as candidates for Clear-Sky samples.
%
%   'tol',t,'minCSfraction',X,'K',k - set by default by SimOptions 'RelTol', 'meteo.minCSfraction',
%       and (indirectly) 'outliers.P'; control convergence tolerance, the threshold to consider
%       if enough clear-samples have been detected, and the threshold for outlier detection.
%
% [C,SP] = FITCLEARSKY(...,'-verbose') - currently does nothing except enable -plot
% [C,SP] = FITCLEARSKY(...,'-plot') - Plot TL and iteration progress (meant for debugging)
%
% See also: PVLMOD_DETECT_CLEAR_TIMES, PVLMOD_CLEARSKY_INEICHEN, LINKETURBIDITY, EFFECTIVESOLARPOSITION

    [opt,varargin] = getflagoptions(varargin,{'-verbose','-plot'});

    opt.nmin = [];
    opt.P = 90;
    opt.minGHI = 20;    % minGHI & minSunEl can & should be much more restrictive than
    opt.minsunel = 8;   % ... SimOptions.minGHI/minSunEl
    opt.limTL = [1.7,6.0];
    
    opt.weights = [0.5 1 1 1 2]; % lookup, CSGHI, CSBNI, GHI, BNI
    
    opt.MaxIter = 10;

    opt.tol = getSimOption('RelTol');
    opt.minCSfraction = getSimOption('meteo.minCSfraction');
    opt.K = (norminv(1-0.5*getSimOption('outliers.P'))/norminv(0.75)-1)/2; % IQR-ranges to reach P
    
    opt = getpairedoptions(varargin,opt,'restchk');

    validateattributes(MD,{'MeteoData'},{'nonempty'});    
    assert(isfield(MD,'GHI'),'FITCLEARSKY currently requires GHI')
    
    if ~isregular(MD.data) || any(~isfield(MD,{'sunel','AMa','ENI'}))
        MD = MeteoData.loadobj(MD); 
    end
    if MD.interval ~= 'i', MD.interval = 'c'; end
    dt = minutes(MD.timestep);
    
    altitude = MD.location.altitude;
    rms = @(x) sqrt(mean(x.^2,'all','omitnan'));

    if isempty(opt.nmin), opt.nmin = round(1440/dt); end  % 24 (weighted) hours
        
    EHI = MD.data.ENI.*max(0,sind(MD.data.sunel));

    B = meteoQC.flagged2nan(MD,'all',{'GHI','BNI'});
    GHI = B.data.GHI;
    GHI(B.dark,:) = NaN;
    GHI(GHI < 0) = NaN;
    if isfield(B,'BNI')
        BNI = B.data.BNI;
        BNI(B.dark,:) = NaN;
        BNI(BNI < 0) = NaN;
    else
        BNI = []; 
    end
    clear B
    
    TL0 = cell(1,3); % {lookup, CSGHI, CSBNI} 
    
    % Start with seed value from lookup tables
    TL0{1} = linketurbidity(MD.location,MD.t);
    % w{1} = opt.weights(1); % w = [0.5 1 1 1 2]; % lookup, CSGHI, CSBNI, GHI, BNI

    % Average-in smoothed-TL calculated from existing CSGHI, CSBNI
    if isfield(MD,'CSGHI') && ~all(contains(getsourceof(MD,'CSGHI'),'ineichen'))
        new = ~contains(getsourceof(MD,'CSGHI'),'ineichen');
        TL0{2} = reverse_clearsky_ineichen(MD.data.CSGHI(:,new),altitude,MD.data.AMa,EHI);
    end
    if isfield(MD,'CSBNI') && ~all(contains(getsourceof(MD,'CSBNI'),'ineichen'))
        new = ~contains(getsourceof(MD,'CSBNI'),'ineichen');
        TL0{3} = BNI_linke_turbidity(MD.data.CSBNI(:,new),altitude,MD.data.AMa,MD.data.ENI);
    end
    
    n = cellfun(@(x) size(x,2),TL0);
    TL0 = [TL0{n > 0}];
    [~,TL0] = minQC(TL0,'sequence',{'num','IQR'},'valid',@(x) x > opt.limTL(1) & x < opt.limTL(2));

    if size(TL0,2) > 1
        w = repelem(opt.weights(1:3),n);
        w = w/sum(w);
        W = sum(w.*isfinite(TL0),2);
        TL0 = sum(TL0.*w,2,'omitnan')./W;
        TL0 = smoothlinketurbidity(TL0,W,MD.data.AMa,round(7*1440/dt),opt.P);  % 7-day average
    end
    
    w = [sum(opt.weights(1:3)),opt.weights(4:5)];
    nGHI = size(GHI,2);
    if ~isempty(BNI)
        nBNI = size(BNI,2);
        TL = [TL0,NaN(MD.Nt,nGHI+nBNI)];
        w = repelem(w,[1,nGHI,nBNI]);
    else
        TL = [TL0,NaN(MD.Nt,nGHI)];
        w = repelem(w(1:2),[1,nGHI]);
    end
    w = w/sum(w);
    
    if isfield(MD,'CSGHI') && (isempty(BNI) || isfield(MD,'CSBNI'))
        CSGHI = mean(MD.data.CSGHI,2);
        if ~isempty(BNI), CSBNI = MD.data.CSBNI; else, CSBNI = []; end
    else
        [CSGHI,CSBNI] = pvlmod_clearsky_ineichen(MD.location,[],TL0,MD.data.sunel,MD.data.AMa,MD.data.ENI);
    end
    N = nnz(CSGHI > 0 & any(GHI > 0,2));
    
    warning_disabler = [];
    if isfield(MD,'clearsky') && nnz(MD.data.clearsky)/N > opt.minCSfraction
        CS = MD.data.clearsky;   
    else
        CS = detect_clear_times(GHI,CSGHI,dt); 
    end
    
    function CS = detect_clear_times(GHI,CSGHI,dt)
        if size(GHI,2) > 1
            [L,H] = bounds(GHI,2,'omitnan');
            GHI = mean(GHI,2,'omitnan');
            fishy = (H-L) > 75; % max. diff in pvlmod_detect_clear_times
        else
            fishy = [];
        end
        
        CS = pvlmod_detect_clear_times(GHI,CSGHI,dt); 
        CS(fishy) = false;
        if isempty(warning_disabler)
            warning_disabler = naptime('detectcleartimes:notoneminute');
        end
    end

    % initialize state for iterative process
    lastCSGHI = CSGHI;
    lastCSBNI = CSBNI;
    best_err = Inf;

    for j = 1:opt.MaxIter
    % Iterative CS search and TL correction based on clear-samples
        
        if j > 1, CS = detect_clear_times(GHI,CSGHI,dt); end
    
        % if j > 1 || nnz(CS)/N < opt.minCSfraction
            TL(:,1) = TL0(:,1);
            TL(:,2:end) = NaN; % reset all but seed

            % Get a smooth Linke-Turbidity series based on GHI(CS) and alpha
            TL(CS,2:1+nGHI) = reverse_clearsky_ineichen(GHI(CS,:),altitude,MD.data.AMa(CS),EHI(CS));

            % Remove outliers based on TL(BNI)
            if ~isempty(BNI)
                TL(CS,2+nGHI:end) = BNI_linke_turbidity(BNI(CS,:),altitude,MD.data.AMa(CS),MD.data.ENI(CS));
            end

            [~,TL] = minQC(TL,'sequence',{'num','IQR'},'valid',@(x) x > opt.limTL(1) & x < opt.limTL(2));
            CS(~any(TL > opt.limTL(1),2)) = false;
            CS(any(TL > opt.limTL(2),2)) = false;
            TL(~CS,:) = NaN;
        % end
        
        if nnz(CS)/N < opt.minCSfraction
            TLs = TL0;
            warning('Failed to detect enough clear-sky samples');
            break;
        end
        
        W = sum(w.*isfinite(TL),2);
        TLs = sum(TL.*w,2,'omitnan')./W;
        TLs = smoothlinketurbidity(TLs,W,MD.data.AMa,opt.nmin,opt.P);

        % Update the CSGHI & CS estimates based on new TL
        [CSGHI,CSBNI] = pvlmod_clearsky_ineichen(MD.location,[],TLs,MD.data.sunel,MD.data.AMa,MD.data.ENI);        
        
        if opt.plot
        % Plot progress, usually for debugging
            if j == 1
                ERR = NaN(opt.MaxIter,4);
                GUIfigure('fitclearsky'); clf();
            end
            subplot(2,1,1); cla(); hold on;
            plot(MD.t,TL0);
            plot(MD.t,TL);
            plot(MD.t,TLs);
            datetick(gca,'x','keeplimits');
            ylabel('Linke Turbidity');
            legend('Seed','Clear-Samples','Smoothed');
            
            ERR(j,1) = rms(lastCSGHI-CSGHI);
            ERR(j,2) = rms(lastCSGHI(CS)-CSGHI(CS));
            ERR(j,3) = rms(CSGHI(CS) - GHI(CS,:));
            ERR(j,4) = rms(TLs(CS)-TL(CS));

            subplot(2,1,2); cla();
            if j > 1
                plot(ERR(1:j,:)./ERR(j,:));
                ylabel('Relative Error');
                xlabel('Iteration');
                legend('CSGHI','CSGHI(CS)','CSGHI-GHI','TLs-TL')
            end
            drawnow()
        end
        
        if isfield(MD,'BNI')
            err = rms([lastCSGHI(CS)-CSGHI(CS);lastCSBNI(CS)-CSBNI(CS)])/...
                  rms([CSGHI(CS) - GHI(CS);CSBNI(CS)-BNI(CS)]);
        else
            err = rms(lastCSGHI(CS)-CSGHI(CS))/rms(CSGHI(CS) - GHI(CS));
        end
        if err < opt.tol || err > best_err
            break; 
        end
        
        lastCSGHI = CSGHI;
        lastCSBNI = CSBNI;
        best_err = min(best_err,err);

        if j == opt.MaxIter
            warning('Failed to converge withing max. allowed iterations');
            break;
        end
    end

    if isfield(MD,'clearsky') && contains(getsourceof(MD,'clearsky'),'detect_clear_times')
        MD = rmfield(MD,'clearsky');
    end
    MD = addsource(MD,'clearsky',CS,'detect_clear_times',true);
    
    if isfield(MD,'TL') && contains(getsourceof(MD,'TL'),{'linketurbidity','fitclearsky'})
        MD = rmfield(MD,'TL');
    end
    MD = addsource(MD,'TL',TLs,'fitclearsky.TL',true);
    MD = getsunpos(MD,0,1);
end

function TLs = smoothlinketurbidity(TL,W0,AMa,NMIN,P)
% Calculate a moving-window weigted average of TL, choosing a window width that ensures the
% equivalent of NMIN full-weight samples for each new averaged value.

    % if nargin < 3, NMIN = 1; end
    % if nargin < 4, P = 50; end
    
    wfun = @(x) exp(-(x-2).^2);  % decreasing weight to points with AM ~= 2
    W = W0.*wfun(AMa);        % an hour with full data at AM < 2 will have W = 1.

    % Get a minimum window-width for areas with low point density.
    N = numel(TL);
    CW = cumsum(W,1,'omitnan');
    pp = mdpwl(-CW,1:N);
    s(:,1) = floor(pp.val(-CW+NMIN/2));
    s(:,2) = ceil(pp.val(-CW-NMIN/2));
    s(s(:,1) < 1,2) = ceil(pp.val(-NMIN));
    s(s(:,2) > N,1) = floor(pp.val(-CW(end)+NMIN));
    s = min(max(1,s),N);
    minwindow = max(diff(s,1,2));

    TLs = movsum(TL.*W,minwindow,'omitnan')./movsum(W,minwindow,'omitnan');
    
    % Offset to median (or other user-selected percentile)
    TLs = TLs - prctile(TLs-TL,P);
end

function TL = BNI_linke_turbidity(BNI,altitude,AMa,ENI)
% Calculate Linke-Turbidity-Factors that would explain Clear-Sky-BNI for a given altitude, air-
% mass, and normal extraterrestrial irradiance, according to Ineichen clear-sky model.
    b = 0.664 + 0.163 ./ exp(altitude.*(-1/8000));
    TL = 1 - log(BNI./(b.*ENI))./(0.09*AMa);
    TL = TL - 0.25*sqrt(max(0,2 - TL));
end

function TL = reverse_clearsky_ineichen(csGHI,altitude,AMa,EHI)
% Calculate Linke-Turbidity-Factors that would explain Clear-Sky-GHI for a given altitude, air-
% mass, and horizontal extraterrestrial irradiance, according to Ineichen clear-sky model.

    fh1 = exp(altitude.*(-1/8000));
    fh2 = exp(altitude.*(-1/1250));
    cg1 = (0.0000509.*altitude+0.868);
    cg2 = 0.0000392.*altitude+0.0387;

    A = cg1.*EHI.*exp(0.01.*(AMa).^(1.8)).*exp(-cg2.*AMa.*(fh1-fh2));
    B = -cg2.*AMa.*fh2;

    %dark = isnan(A) | isnan(B) | A <= 0;
    %A(dark) = 0; B(dark) = -Inf;
    
    TL = (log(csGHI)-log(A))./B;
    TL(csGHI < 0 | A < 0) = NaN;
end
