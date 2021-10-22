function [Kd,Kn,mdl] = diffuse_fraction(varargin)
% [KD,KN] = DIFFUSE_FRACTION(KT,...,MODEL) - Wrapper for separation models
%
%   'orgill' - Orgill-Hollands (1977) model [1]. See PVL_ORGILL-HOLLANDS.
%   'erbs' - Erbs et al. (1982) model [2]. See PVL_ERBS.
%   'reindl1' - Reindl et al. (1990) 1-parameter model [3]. See PVL_REINDL_1.
%   'reindl2' - Reindl et al. (1982) 2-parameter model [3]. See PVL_REINDL_2.
%   'disc' - Maxwell (1987) DISC model [4]. See PVLMOD_DISC.
%   'dirint' - Perez et al. (1992) DIRINT model [5]. See PVLMOD_DIRINT.
%   'engerer2' - Engerer2
%
% Inputs: Optional positional arguments after Kt are read as: { SZA, Patm, TPW, timeseries, Ktc, 
%   hourangle, AM }. The same can be provided as name,value pairs, e.g. 'SZA',Z
%
%   Kt - array of clearness indices (Kt = GHI / Extraterrestrial-Horizontal-Irradiance)
%   SZA - scalar or size(Kt) array of true effective solar zenith angles [degrees], required for
%       reindl2, DISC and DIRINT models.
%   Patm - scalar or size(Kt) array of atmospheric pressure values [Pa], default 101325 Pa,
%       only required for DISC and DIRINT models.
%   TPW - scalar or size(Kt) array of precipitable water content [mm]. Only used by DIRINT model.
%   timeseries - boolean is-time-series flag in PVLMOD_DIRINT
%   Ktc - clear-sky clearness index (CSGHI/EHI) for engerer2 and sot3 models.
%   hourangle - solar hour angle (degrees) used to get AST in Engerer2model
%   AM - Absolute air mass, not required if SZA and Patm are also provided.
%
%   MODEL - character key for model selection (see above). The default is chosen based on the
%       available arguments, with priority for the last models of the list above.
%
% Output:   
%   Kd - array of diffuse fraction values (Kd = DHI / GHI)
%   Kn - array of direct clearness indices (Kn = DNI / GHI)
%
% Sources: 
%
%   [1] Orgill JF, Hollands KGT. Correlation equation for hourly diffuse ratiation on a horizontal
%       surface. Solar Energy 1977;19:357-9 
%   [2] D. G. Erbs, S. A. Klein and J. A. Duffie, Estimation of the diffuse radiation fraction for
%       hourly, daily and monthly-average global radiation, Solar Energy 28(4), pp 293-302, 1982.
%   [3] Reindl DT, Beckman WA, Duffie JA. Diffuse fraction correlations. Solar Energy 1990;45:1-7.
%   [4] Maxwell, E. L., "A Quasi-Physical Model for Converting Hourly Global Horizontal to Direct 
%       Normal Insolation", Technical Report No. SERI/TR-215-3087, Golden, CO: Solar Energy  
%       Research Institute, 1987.
%   [5] Perez, R., P. Ineichen, E. Maxwell, R. Seals and A. Zelenka, (1992). Dynamic Global-to-
%       Direct Irradiance Conversion Models.  ASHRAE Transactions-Research Series, pp. 354-369
%   [6] Code for Engerer2 was adapted from Yang & Boland, 2019. “Satellite-Augmented Diffuse Solar 
%       Radiation Separation Models.” Journal of Renewable and Sustainable Energy 11 (2): 023705.
%
% TODO: translate other models from Yang & Boland. 2019, accept 'ensemble' as model?
%
% See also PVLMOD_DIRINT, PVLMOD_DISC, PVL_ORGILL_HOLLANDS PVL_REINDL_1 PVL_REINDL_2 PVL_ERBS 

    narginchk(1,Inf);
    argnames = {'kt','sza','patm','tpw','timeseries','ktc','hourangle','am','dt'};
    
    % last models will be picked first, if all arguments are available
    MDL = { 'erbs'    ,{'kt'}, 'kd', @erbs;
            'orgill'  ,{'kt'}, 'kd', @orgill_hollands;
            'reindl1' ,{'kt'}, 'kd', @reindl1;
            'reindl2' ,{'kt','sza'}, 'kd', @reindl2;
            'disc'    ,{'kt','am'}, 'kn', @pvlmod_disc;
            'sot2'    ,{'kt','sza'}, 'kd', @skartveit_olseth_tuft;
            'sot3'    ,{'kt','sza','ktc'}, 'kd', @skartveit_olseth_tuft;
            'sot4'    ,{'kt','sza','ktc','dt'}, 'kd', @skartveit_olseth_tuft;
            'dirint'  ,{'kt','sza','patm','tpw','timeseries'}, 'kn', @pvl_dirint;
            'engerer2',{'kt','sza','hourangle','ktc'}, 'kd', @engerer2};
        
    MDL = cell2struct(MDL',{'name','input','output','fcn'});
    
    if ischar(varargin{end})
    % Explicit model request
        [~,ic] = parselist(varargin{end},{MDL.name},'model');
        model = MDL(ic);
        varargin(end) = [];
    else
        model = [];
    end
    
    % parse both positional & name-value pair arguments
    args = getpairedoptions(varargin,argnames,'dealrest');
    emptyfields = fieldnames(args);
    emptyfields(~cellfun(@(f) isempty(args.(f)),emptyfields)) = [];
    args = rmfield(args,emptyfields);
    
    parsestruct(args,{'kt'},'opt',{'sza','patm','tpw','ktc','hourangle','am'},'-n','-r','-c');
    if isfield(args,'timeseries')
        validateattributes(args.timeseries,{'numeric','logical'},{'binary','scalar'});
    end  
    if isfield(args,'dt')
        validateattributes(args.dt,{'numeric','duration'},{'real','scalar','positive','nonzero'});
        if isduration(args.dt), args.dt = minutes(args.dt); end
    end
    
    used = {}; % args used indirectly will be placed here, to skip 'ignored argument' warning
    if ~isfield(args,'am') && all(isfield(args,{'sza','patm'}))
        args.am = pvl_relativeairmass(args.sza,'kasten1966').* args.patm./101325;
        used = {'sza','patm','am'}; 
    end
    
    if isempty(model)
    % pick the last (~best) model for which required arguments are available
        couldwork = arrayfun(@(mdl) all(ismember(mdl.input,fields(args))),MDL);
        model = MDL(find(couldwork,1,'last'));
        
        if nargout < 3
            warning('diffuse_fraction:defaultmdl','Using default diffuse model: %s',model.name);
        end
    end
    mdl = model.name;
    
    if ~isfield(args,'patm') && ismember('patm',model.input)
        args.patm = 101325;
        warning('Assuming sea-level pressure');
    end
    
    Kt = args.kt;
    available = fieldnames(args);
    args = struct2cell(args);

    [reallythere,idx] = ismember(model.input,available);
    assert(all(reallythere),shortliststr(model.input(~reallythere),'Missing model argument'));

    ignored = setdiff(available,[used,available(idx)']);
    if ~isempty(ignored)
        warning('diffuse_fraction:ignored','''%s'' model ignores %s',...
            model.name,shortliststr(ignored));
    end
    
    args = args(idx(reallythere));
    kx = model.fcn(args{:});
    switch model.output
        case 'kd', Kd = kx; if nargout > 1, Kn = Kd2Kn(Kd,Kt); end
        case 'kn', Kn = kx; Kd = Kn2Kd(Kn,Kt);
        otherwise
            error('you really should not be here');
    end
end

function Kd = Kn2Kd(Kn,Kt)
    Kd = 1 - Kn./Kt;
    Kd(Kn <= 0 | Kt <= 0) = 1;
end

function Kn = Kd2Kn(Kd,Kt)
    Kn = Kt.*(1-Kd);
    Kn(Kt <= 0 | Kd >= 1) = 0;
end

function kd = reindl1(kt)
    kt(kt > 1) = 1;
    kt(kt < 0) = 0;
    kd = interp1([0,0.3,0.78,1],[1,0.9456,0.147,0.147],kt);
end

function kd = reindl2(kt,Z)

    [kt,Z] = compatiblesize(kt,Z);
    kd = NaN(size(kt));

    i = kt <= 0.30;
    kd(i) = 1.02 - 0.254*kt(i) + 0.0123*cosd(Z(i));
        
    i = kt > 0.30 & kt < 0.78;
    kd(i) = 1.4 - 1.749*kt(i) + 0.177*cosd(Z(i));
    kd(i) = min(max(kd(i),0.1),0.97);
        
    i = kt >= 0.78;
    kd(i) = 0.486*kt(i) - 0.182*cosd(Z(i));
end

function kd = orgill_hollands(kt)

    kd = NaN(size(kt));

    i = kt <= 0.35;
    kd(i) = 1.0 - 0.249*kt(i);
        
    i = kt > 0.35 & kt < 0.75;
    kd(i) = 1.577 - 1.84*kt(i);
        
    i = kt >= 0.75;
    kd(i) = 0.177;
end

function kd = erbs(kt)

    kd = NaN(size(kt));

    i = kt <= 0.22;
    kd(i) = 1 - 0.09.*kt(i);

    i = kt>0.22 & kt<=0.8;     
    kd(i) = 0.9511 +kt(i).*(-0.1604 + kt(i).*(4.388 +kt(i).*(-16.638 + kt(i).*12.336)));

    i = kt > 0.8;
    kd(i) = 0.165;
end

function kd = engerer2(Kt, Z, w, Ktc)

    AST = w/15 + 12;
    delta_Ktc = Ktc - Kt;
    Kde = max(0,1 - Ktc./Kt);

    C = 4.2336e-2;
    beta0 = -3.7912; beta1 = 7.5479; beta2 = -1.0036e-2;
    beta3 = 3.1480e-3; beta4 = -5.3146; beta5 = 1.7073;

    kd = C + (1-C)./(1+exp(beta0+beta1*Kt + beta2*AST + beta3*Z + beta4*delta_Ktc)) + beta5*Kde;
end

function d = skartveit_olseth_tuft(k,Z,Ktc,dt)

    if nargin < 3, Ktc = []; end
    if nargin < 4, dt = []; end
    if ~isempty(dt)
       if mod(60,dt) ~= 0 || dt <= 60
          warning('Cannot use timesteps that are not an integer divisor of 1 hour');
          dt = [];
       end
    end

    if ~isempty(Ktc)
        [k,Z,Ktc] = compatiblesize(k,Z,Ktc);
        r = k/Ktc;
        
        if ~isempty(dt) && (numel(Z) < 180/dt || ~isvector(Z))
            warning('Using average variability correction');
            dt = [];
        end

        s3 = variability_index(r,dt);
    else
       [k,Z] = compatiblesize(k,Z);
    end

    % Invariable hours
    h = 90-Z;
    [k,h] = compatiblesize(k,h);
    
    d = ones(size(k)); % d = 1 for k < 0.22, eq. (4)

    d1 = 0.07+0.046*(90-h)./(h+3); % (6d)
    d1(h < 1.4) = 1;
    k1 = 0.83 - 0.56*exp(-0.06*h); % (6b)
    k2 = 0.95*k1; % (6c)
    K = 0.5*(1+sin(pi*(k2-0.22)./(k1-0.22)-pi/2)); % (6a)
    d2 = 1 -(1-d1).*(0.11*sqrt(K)+0.15*K+0.74*K.^2); % (5)

    alpha = cscd(h).^0.6; % (9b)
    kbmax = 0.81.^alpha;  % (9a)
    kmax = (kbmax + d2.*k2./(1-k2))./(1+d2.*k2./(1-k2));
    
    % For 0.22 < k < k2
    ib = (k > 0.22) & (k <= k2);
    K(ib) = 0.5*(1+sin(pi*(k(ib)-0.22)./(k1(ib)-0.22)-pi/2)); % (6a)
    d(ib) = 1 -(1-d1(ib)).*(0.11*sqrt(K(ib))+0.15*K(ib)+0.74*K(ib).^2); % (5)
    
    % For k2 < k < kmax
    ib = (k > k2) & (k <= kmax);
    d(ib) = d2(ib).*k2(ib).*(1-k(ib))./(k(ib).*(1-k2(ib))); % (7)

    % For k > kmax
    ib = (k > kmax);
    dmax = d2(ib).*k2(ib).*(1-kmax(ib))./(kmax(ib).*(1-k2(ib))); % (11)
    d(ib) = 1 - kmax(ib).*(1-dmax)./k(ib); % (10)

    if ~isempty(Ktc)
    % Variability

        kx = 0.56-0.32*exp(-0.06*h); % (13a)

        D = zeros(size(k));  % for k < 0.14 or k > kx + 0.71 (12c)
        ib = (k > 0.14) & k <= kx;
        kL = (k(ib) - 0.14)./(kx(ib) - 0.14);  % (13b)
        D(ib) = -3*kL.^2 * (1-kL) * s3(ib).^1.3;  % (12a)

        ib = (k > kx) & k <= kx + 0.71;
        kR = (k(ib)-kx(ib))/0.71;  % (13c)
        D(ib) = 3*kR * (1-kR).^2 * s3(ib).^0.6;  % (12b)

        d = d + D;
    end
end

function s3 = variability_index(r,dt)
% Hourly variability index according to Skartveit, Olseth & Tuft, 1998.

    MINAVAIL = 0.5; 
    
    s3 = NaN(size(r));
   
    if isempty(dt)
    % use average variability index (eqns. 3a-3b)
        ib = (r <= 1.014);
        s3(ib) = 0.021 + 0.387*r(ib) - 0.231*r(ib).^2 - 0.13./exp(((r(ib)-0.931)/0.134).^(2*0.834));
        s3(~ib) = 0.12 + 0.65*(r(~ib)-1.04);
    else
        % There is no guideline in Skartveit et al. on what to do for timeseries of higher
        % resolution, but a one-hour moving window average should be equivalent.
        
        n = 60/dt;
        R = movmean(r(:),n,'omitnan');
        R(movmean(isfinite(r(:)),n) < MINAVAIL) = NaN;
        R = [circshift(R,n),R,circshift(R,n)];
        R(1:n,1) = NaN;
        R(end-n+1:end,3) = NaN;
        ib = all(isfinite(R),2);
        s3(ib) = hypot(R(ib,2) - R(ib,1),R(ib,2) - R(ib,3))/sqrt(2);

        % Use the absolute difference if only two samples are available
        ib = all(isfinite(R(:,1:2),2));
        s3(ib) = abs(R(:,1)-R(:,2));
        ib = all(isfinite(R(:,2:3),2));
        s3(ib) = abs(R(:,2)-R(:,3));
        
        ib = ~isfinite(s3) & isfinite(r);
        s3(ib) = variability_index(r(ib));
    end
end
