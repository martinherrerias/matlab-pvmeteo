function [CSGHI,CSBNI,CSDHI] = clearsky_ESRA(Loc,t,TL,sunel,mdl)
% [CSGHI,CSBNI,CSDHI] = clearsky_ESRA(Loc,t,[TL,sunel,mdl]) - Modified ESRA clear-sky model,
%   according to Rigollier et al. 2000 and Geiger et al. 2002.
%
% Input Parameters:
%   LOC - a struct with scalar fields 'latitude', 'longitude' (degrees N, E), and 'altitude' (MASL)
%   T - [Nt,1] vector of UTC DATENUM time-steps-
%   TL - (optional) [Nt,1] vector of (AM2) Linke-Turbidity values (IMPORTANT! see note 0)
%   SUNEL - (optional) [Nt,1] apparent SOLAR-elevation angle (degrees)
%   MDL - 'Rigollier' uses the site elevation correction in [1], for compatibility with [5]. By
%       default ('Geiger'), the correction of [2] is used, since they claim, "The formula for
%       the Rayleigh optical thickness in Rigollier et al. behaves incorrectly with terrain 
%       altitude in the original model."
%
% Output:   
%   CSGHI - the modeled global horizonal irradiance in W/m²
%   CSBNI - the modeled direct normal irradiance in W/m²
%   CSDHI - the modeled diffuse horizonal irradiance in W/m²
%
% Notes:
%   (0) To be consistent with PVL_CLEARSKY_INEICHEN, Linke Turbidity values taken as input
%       are assumed to have a factor 0.8662, i.e. TL = 0.8662 TL* where TL* is defined in
%       agreement with [1]-[2].
%
%   (1) Both the code from Antonazas-Torres et al. (2019) [7] as well as the C code in [6]
%       use inverse-square-weighting instead of linear-interpolation for site-elevation correction.
%       Based on personal communication with L. Wald [3], linear interpolation should be used, as
%       documented in [2].
%
%   (2) A. Equation 6.a of [2] reads: 1/dR = corr_dR / (6.62 + ... ). This seems to be a typo,
%          as the code in [6] uses corr_dR in the denominator.
%       B. Both equation 6.b in [2] as the the code in [6] generate an abrupt discontinuity
%          at m = 20 when p/p0 < 1. It seems to arise from the use of "m" for both relative- or 
%          absolute-air-mass under different circumstances. 
%       The issues 2.A. and 2.B. were corrected by using:
%
%       1/dR = corr_dR x (6.62 + 1.93 m + ... ) for m <= 20
%       1/dR = ( 10.4 + 0.718·m·(p/p0) )/(p/p0) for m > 20
%       TRB = exp(-0.8662·TL·m·dR)
%
%       For m = relative-air-mass. Note that there is still a small discontinuity at m = 20
%       (for all values of p/p0) but it has a minor effect.
%
%   (3) Equations 9 and 11 in [1] show polynomials in TL, but the code [6] uses pressure-
%       corrected TLAM2_corr = TL·p/p0. This function uses the latter.
%
% References:
% [1] C. Rigollier, O. Bauer, and L. Wald, “On the clear sky model of the ESRA — European Solar 
%   Radiation Atlas — with respect to the heliosat method,” Solar Energy, vol. 68, no. 1, 
%   pp. 33–48, Jan. 2000, doi: 10.1016/S0038-092X(99)00055-9.
% [2] M. Geiger, L. Diabaté, L. Ménard, and L. Wald, “A web service for controlling the quality
%   of measurements of global solar irradiation,” Solar Energy, vol. 73, no. 6, pp. 475–480, 
%   Dec. 2002, doi: 10.1016/S0038-092X(02)00121-4.
% [3] L. Wald, Personal Communication. Sep. 20, 2019.
% [4] J. Remund, L. Wald, M. Lefèvre, T. Ranchin, and J. Page, “Worldwide Linke turbidity 
%   information,” ISES Solar World Congress, p. 14.
% [5] https://github.com/JamieMBright/clear-sky-models
% [6] Clear Sky library: Software library for implementing the ESRA clear sky model providing 
%   estimates of the surface solar irradiance under clear sky
%   https://www.oie.minesparis.psl.eu/Valorisation/Outils/Clear-Sky-Library/
%   https://www.oie.minesparis.psl.eu/Donnees/data12/1205-OIE_csmodels_lib.c
% [7] F. Antonanzas-Torres, R. Urraca, J. Polo, O. Perpiñán-Lamigueiro, and R. Escobar, 
%   “Clear sky solar irradiance models: A review of seventy models,” Renewable and Sustainable
%   Energy Reviews, vol. 107, pp. 374–387, Jun. 2019, doi: 10.1016/j.rser.2019.02.032.
%
% See also: LINKETURBIDITY, PVL_CLEARSKY_INEICHEN, PVL_CLEARSKY_HAURWITZ, EFFECTIVESOLARPOSITION

    if nargin == 0, test(); return; end

    narginchk(2,5);
    parselocation(Loc,'optional',{'TimeZone','name'});
    t = parsetime(t);
    t = t(:);
       
    if nargin < 3, TL = []; end
    if nargin < 4, sunel = []; end
    if nargin < 5 || isempty(mdl), mdl = 'geiger'; end
    
    if ~isempty(t)
    % Get any/all missing time-dependent variables
        if isempty(sunel)
            [~,~,sunel] = pvlmod_ephemeris(t,Loc);
            sunel = reshape(sunel,size(t));
        end
        if isempty(TL), TL = linketurbidity(Loc,t); end   
    else
        narginchk(4,5);
    end
    validateattributes(sunel,{'numeric'},{'real','vector','size',[NaN,1]});
    validateattributes(TL,{'numeric'},{'real','vector','size',[NaN,1]});
    compatiblesize(t,sunel,TL);
    
    TL = TL/0.8662; % see note (0)
    
    cosSZA = sind(sunel);
    
    % m = pvl_relativeairmass(90-sunel,'kastenyoung1989');
    m = 1./(cosSZA+0.50572.*max(0,6.07995 + sunel).^-1.6364); % clip only at -6°
    m(sunel < 0) = NaN;
    
    % AMa = exp(-Loc.altitude/8434.5).*AMr;
    ENI = solarposition.extraradiation(t,1367);

    p_p0 = exp(-Loc.altitude/8434.5);
    
    switch lower(mdl)
    case {'kasten','rigollier','2000'}
    
        m = m.*p_p0;
        
        C = [6.6296,1.7513,-0.1202,0.0065,-0.00013];
        dR = 1./polyval(fliplr(C),m);
        dR(m > 20) = 1./( 10.4 + 0.718*m(m > 20) );
        
    case {'page','geiger','2002'}

        % See note (1) on interpolation
        corr_C = [1.682190, -0.030590, 0.00089;  % p/p0 = 0.5
                  1.248174, -0.011997, 0.00037;  % p/p0 = 0.75
                  1,0,0];                        % p/p0 = 1.0

        W = interpmatrix(p_p0,[0.5,0.75,1]','extrap','nearest');
        corr_C = full(W)*corr_C;
        corr_dR = polyval(fliplr(corr_C),m);

        C = [6.625928, 1.92969,-0.170073,0.011517,-0.000285];
        dR = (1./polyval(fliplr(C),m))./corr_dR;                % see note (2A)
        dR(m > 20) = p_p0./( 10.4 + 0.718*m(m > 20)*p_p0 );     % see note (2B)
    end

    % direct beam irradiance
    CSBNI = max(0,ENI.*exp(-0.8662*TL.*m.*dR));
    CSBNI(sunel < 0) = 0;
    
    TL_p2 = (TL.*p_p0).^[0,1,2]; % see note (3)
   
    % the diffuse transmission function at zenith
    TRD = TL_p2*[-1.5843e-2,3.0543e-2,3.797e-4]';
    
    % Set of quadratic functions [A0,A1,A2]
    A = [0.26463,-6.1581e-2, 3.1408e-3;  % a0
         2.0402 , 1.8945e-2,-1.1161e-2;  % a1
         -1.3025, 3.9231e-2, 8.5079e-3]; % a2
    A = (TL_p2)*A';
    f = TRD.*A(:,1) < 2.1e-3; % A0 = 2.1e-3/Td for Td·A0 < 2.1e-3
    A(f,1) = 2.1e-3./TRD(f);
    FD = sum(A.*cosSZA.^[0,1,2],2); % diffuse angular function

    % diffuse horizontal irradiance
    CSDHI = ENI.*max(0,TRD.*FD);
    
    % global horizontal irradiance
    CSGHI = CSBNI.*cosSZA+CSDHI;
end

function test()
    
    N = 1;
    
    lat = (rand(N,1)*2-1)*70;
    lon = (rand(N,1)*2-1)*180;
    alt = rand(N,1)*3000;
    
    t0 = datenum(2019,1,1)+(0.5:1440)'/1440 + solarposition.mdoy(1:2:12)-1;
    t0 = datetime(t0(:),'convertfrom','datenum');

    for j = 1:N
        Loc = struct('latitude',lat(j),'longitude',lon(j),'altitude',alt(j));
        t = t0; t.TimeZone = checktimezone(round(lon(j)/15));
        TL = linketurbidity(Loc,t);
        [~,sunel] = pvlmod_ephemeris(t,Loc);
        [CSGHI,CSBNI,~] = clearsky_ESRA(Loc,t,TL,sunel,'geiger');
        [CSGHI2,CSBNI2,~] = pvlmod_clearsky_ineichen(Loc,t,TL,sunel,[],[],false);
        
        GUIfigure('clearsky_ESRA_test'); clf(); hold on;
        title(strjoin({deg2dms(lat(j),'NS','precision','1d'),...
                       deg2dms(lon(j),'EW','precision','1d'),...
                       num2str(alt(j),'%0.0f mAMSL')},', '));
        plot(CSGHI);
        plot(CSBNI);
        % plot(CSDHI);
        plot(CSGHI2,'--');
        plot(CSBNI2,'--');
        % plot(CSDHI2,'--');
        legend('CSGHI_{ESRA}','CSBNI_{ESRA}','CSGHI_{Ineichen}','CSBNI_{Ineichen}');
        
        if j < N, pause(); end
    end
end