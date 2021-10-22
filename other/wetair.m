classdef wetair
% Set of APPROXIMATE Arden Buck (1981) equations for Air-Water Psychrometrics. 
% Unless otherwise stated, the provided approximations and ideal gas relations should ensure 
% an error below 2% for air between -20°C and 50°C at atmospheric pressure. See the original
% publication [1] for more precise approximations.
%
% [1] A. L. Buck, “New Equations for Computing Vapor Pressure and Enhancement Factor”,
%   J. Appl. Meteor., vol. 20, no. 12, pp. 1527–1532, Dec. 1981.
%
% (Ɔ) Martín Herrerías Azcué, 2009-2019.

methods (Static)

    function p = alt2pres(altitude)
    % P = ALT2PRES(MASL) - return approx. pressure (Pa) from altitude (MASL)
    %   Copied from PVLIB (PVPVMC), for convenience.
    %
    % Assumptions include:
    %   Base pressure = 101325 Pa
    %   Temperature at zero altitude = 288.15 K
    %   Gravitational acceleration = 9.80665 m/s²
    %   Lapse rate = -6.5E-3 K/m
    %   Gas constant for air = 287.053 J/(kg·K)
    %   Relative Humidity = 0%
    %
    % References:
    %   "A Quick Derivation relating altitude to air pressure" from Portland
    %   State Aerospace Society, Version 1.03, 12/22/2004.

        p = 100* ((44331.514 - altitude)/11880.516).^(1/0.1902632);
    end

    function pw = vpres(Tdb, RH, Patm)
    % PW = VPRES(TDB,RH,PATM) - returns the partial vapor pressure [Pa] of moist air at temperature
    %   Tdb [°C], relative humidity RH [0,1], and (optionally) atmospheric pressure PATM [Pa].
    %   Based on equations ew4 and fw4 of Arden-Buck (1981) for moist air (on water). Ensures
    %   error < 0.23% between -40°C and 50°C
    
        if nargin < 3 || isempty(Patm), Patm = 101325; end
    
        f = 1.00072 + (Patm/100).*(3.20e-6 + 5.9e-10*Tdb.^2); % enhancement factor (eq. fw4 )
        es = 6.1121 * exp((18.729 - Tdb/227.3) .* Tdb ./ (Tdb + 257.87)); % (eq. ew4)
        
        pw = 100 * f .* es .* RH;
    end

    function DP = dewpoint(pw)
    % DP = DEWPOINT(PW) - returns dew-point temperature [°C] of air with a water vapor pressure
    %   PW [Pa] according to equation w1 of Arden-Buck (1981) for air between -20°C and 50°C (<1% Error)
    
        DP = 240.97 * log(pw / 613.65) ./ (17.502 - log(pw / 613.65));
    end

    function w = mixratio(pw, Patm)
    % W = MIXRATIO(PW,PATM) - returns the mix-ratio (specific humidity) [kg/kg] of wet air at
    %   pressure PATM and water vapor pressure PW.
    
        if nargin < 2, Patm = 101325; end
        
        w = 0.62198 * pw ./ (Patm - pw);
    end
    
    function w = mixratio2(Tdb, Twb, Patm)
    % W = MIXRATIO2(TDB,TWB,PATM) - returns the mix-ratio (specific humidity) [kg/kg] of wet air at
    %   dry- and wet-bulb temperatures TDB [°C] and TWB[°C] respectively, at pressure PATM.
        
       if nargin < 3, Patm = 101325; end

       pw = wetair.vpres(Twb,1);
       w_wb = 0.62198 * pw ./ (Patm - pw); % mix-ratio of saturated air at Twb

        % From wet-bulb definition: h(Tdb,w) = h(Twb,xwb)-(xwb-w)·Cp·Twb, from which:
        w = (1.006 * (Twb - Tdb) + w_wb .* (2501 - (4.18 - 1.84) * Twb)) ./ ...
            (1.84 * Tdb - 4.18 * Twb + 2501);
    end
    
    function h = enthalpy(Tdb, w)
    % H = ENTHALPY(TDB,W) - return the enthalpy of wet air [kJ/kg] at dry-bulb temperature TDB[°C],
    %   and with a mix-ratio of W [kg/kg].
    
        h = 1.006 * Tdb + w .* (1.84 * Tdb + 2501);
    end

    function rho = density(Tdb,RH,Patm)
    % RHO = DENSITY(TDB,RH,PATM) - return the density of wet-air [kg/m³] at dry-bulb temperature
    %   TDB[°C], relative humidity RH, and pressure PATM [Pa].
    
        if nargin < 3, Patm = 101325; end

        pw = wetair.vpres(Tdb,RH);
        % rho = (1/286.9) * (Patm - pw)./(Tdb + 273.15); ??
        rho = (1/287.05) * (Patm - 0.378*pw)./(Tdb + 273.15);
    end
    
     function HR = relhum(Tdb, w, Patm)
    % RH = RELHUM(TDB,W,PATM) - return the relative-humidity of air with temperature TDB[°C] and
    %   mix-ratio W [kg/kg], at pressure PATM.
    
        if nargin < 3, Patm = 101325; end

        pw = Patm .* w ./ (w + 0.62198);        % vapor pressure
        HR = pw ./ wetair.vpres(Tdb, 1);        % relative humidity
    end

    function HR = relhum2(Tdb, Twb, Patm)
    % RH = RELHUM2(TDB,TWB,PATM) - return the relative-humidity of air with dry- and wet-bulb 
    %   temperatures TDB[°C] and TWB[°C], respectively, at pressure PATM.
    
        if nargin < 3, Patm = 101325; end

        w = wetair.mixratio2(Tdb, Twb, Patm);
        pw = Patm .* w ./ (w + 0.62198);        % vapor pressure
        HR = pw ./ wetair.vpres(Tdb, 1);        % relative humidity
    end

    function Twb = wetbulb(Tdb, RH, Patm)
    % TWB = WETBULB(TDB,RH,PATM) - return the wet-bulb temperature [°C] of air at dry-bulb 
    %   temperature TDB [°C] and relative humidity RH [0,1] at pressure PATM [Pa].

        IT_MAX = 100;
        ERR_MAX = 0.01;
        
        if nargin < 3, Patm = 101325; end

        % enthalpy of air at known Tdb and RH:
        pw = wetair.vpres(Tdb, RH);
        w = wetair.mixratio(pw, Patm);
        h = enthalpy(Tdb, w);

        Twb = Tdb .* RH;   % Seed value of Twb for iteration
        
        for i = 1:IT_MAX
            % get enthalpy h(Tdb,Twb) with estimated Twb
            x2 = wetair.mixratio2(Tdb, Twb, Patm); 
            h2 = enthalpy(Tdb, x2);

            % fix Twb based on error in enthalpy
            Twb2 = (h ./ h2) .* Twb;         
            err = abs(Twb2 - Twb);
            if err <= ERR_MAX, break; end
            Twb = Twb2;
        end
    end

    function tpw = precipitablewater(Tdb,RH,varargin)
    % PPW = PRECIPITABLEWATER(TDB,RH,[PATM]) - estimate total atmospheric precipitable water [mm]
    %   or [kg/m²] based on temperature TDB [°C], relative humidity [0,1] and pressure PATM [Pa].
    %   
    % [1] Gueymard, C., 1994. Analysis of monthly average atmospheric precipitable water and 
    % turbidity in Canada and Northern United States. Solar Energy 53, 57–71.

        T = Tdb + 273.15;
        theta = T./273.15;
        
        % Originally uses Gueymard (1993) approximation:
        % es = exp(22.330 - 49.140.*(100./T) - 10.922.*(100./T).^2 - 0.39015.*T./100); % mbar
        
        es = wetair.vpres(Tdb,1,varargin{:}); % Pa
        rho_v = 2.16675*es.*RH./T; % water vapor density [g/m^3], (2.16 = 1000/R)
        
        % Apparent water vapor scale height (km) [1]
        H_v = 0.4976 + 1.5265.*theta + exp(13.6897.*theta - 14.9188.*theta.^3);
        
        tpw = H_v.*rho_v; % ppw-mm (kg/m²)
    end
end
end


