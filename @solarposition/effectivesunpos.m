function [SP,C] = effectivesunpos(Loc,t,varargin)
% [SP,C] = EFFECTIVESUNPOS(LOC,T,[PATM,TA,TL]) - Calculate the effective, apparent
%   solar position for finite intervals [T - dT/2, T + dT/2] at location LOC.
%
%   Following the recommendation of Blanc & Wald (2016), the effective solar zenith Z is such
%   that it supports the closure of a clear-sky model: cos(Z) = (CSGHI - CSDHI)/CSBNI. 
%   In this case, PVLMOD_CLEARSKY_INEICHEN is used, with optional input Linke-Turbidity TL.
%
%   Solar position is calculated by PVLMOD_SPA (considering optional ambient temperature TA and 
%   pressure PATM for refraction). The clear-sky model is evaluated at 1-minute resolution, and
%   'integrated' (downsampled) to the time-series resolution. The solar position for daylight
%   intervals is then updated to match the point for which cos(Z) = (CSGHI - CSDHI)/CSBNI. Note
%   that clear-sky irradiance values are a subproduct, returned in structure CS.
%
% INPUT:
%   LOC - location structure {latitude,longitude,altitude} as required by PVLMOD_SPA
%   T - vector of DATENUM timesteps, expected to be (quasi) regular
%   PATM, TA - (optional) vectors of atmospheric pressure (Pa) and temperature (°C),
%       to account for refraction near the horizon. Defaults are 1.01e5 Pa and 10°C.
%   TL - (optional) vector of Linke-Turbidity values, for PVLMOD_CLEARSKY_INEICHEN
%
%   ..,'step',dt, - start with known time-step dt (see PARSETIME)
%   ..,'interval','x' - use summarization other than 'c' (default) see CHECKSUMMARIZATION.
%
% OUTPUT: 
%   SP - structure with fields {Az,El,w,dec} - vectors of azimuth (N2E) apparent-effective 
%       elevation, hour angle, and declination (all in degrees)
%   C - structure with fields {ENI,CSGHI,CSBNI,CSDHI} - vectors of calculated extraterrestrial 
%       normal irradiance, and clear-sky-global, -beam, and -diffuse irradiances, respectively.
%
% NOTES:
%   Notice that for time-steps < 1.5 min, EFFECTIVESUNPOS directly returns the solar 
%   position at the center of each interval. 
%   
%   For performance, PVLMOD_SPA is calculated at the original resolution, and the 1-minute solar
%   positions required for the clear-sky model are interpolated on hour-angle and declination.
%   This causes errors of O(1e-5) in cos(Z) for hourly intervals & lower for shorter intervals.
%
%   An alternative estimate is SOLARPOSITION.MEANPOSITION, which uses the analytical integral of
%   the top-of-atmosphere horizontal irradiance, cos(Z) = integral(cos(z)·ENI,..)/ENI.
%
% REFERENCES:
%   Blanc, P., Wald, L., 2016. On the effective solar zenith and azimuth angles to use 
%   with measurements of hourly irradiation. Advances in Science and Research 13, 1–6.
%   https://doi.org/10.5194/asr-13-1-2016
%
% See also: PVLMOD_SPA, SOLARPOSITION.MEANPOSITION, SOLARPOSITION.FIXAZIMUTH, ...
%           LINKETURBIDITY, PVLMOD_CLEARSKY_INEICHEN
        
    narginchk(2,14);
    Loc = parselocation(Loc,'-soft');

    opt.step = [];
    opt.interval = 'c';
    opt.Patm = pvl_alt2pres(Loc.altitude);
    opt.Ta = 10;
    opt.TL = [];
    [opt,varargin] = getpairedoptions(varargin,opt);
    assert(numel(varargin) <= 4,'Unrecognized arguments');
    if ~isempty(varargin), Patm = varargin{1}; else, Patm = opt.Patm; end
    if numel(varargin) > 1 && isempty(varargin{2}), Ta = varargin{2}; else, Ta = opt.Ta; end
    if numel(varargin) > 2 && isempty(varargin{3}), TL = varargin{3}; else, TL = opt.TL; end
    
    % Bring labels to center of interval
    opt.interval = checksummarization(opt.interval);
    if opt.interval == 'i'
        t = parsetime(t,'step',opt.step,'interval','i');
        dt = 0;
    else
        [t,dt] = parsetime(t,'step',opt.step,'interval',{opt.interval,'c'});
    end
    
    [t,Patm,Ta] = compatiblesize(t,Patm,Ta);
    if ~isempty(TL), [~,TL] = compatiblesize(t,TL); end

    % Start with refraction-corrected (apparent) solar position at the center of interval
    [SP.Az,~,SP.El,SP.w,SP.dec,r] = pvlmod_spa(t,Loc,Patm,Ta);
    
    C.ENI = single(1361.6./r.^2); % Extraterrestrial Normal Irradiance
        
    k = round(minutes(dt)); 
    if k <= 1
    % For dt >= 1 min, just use the center of the interval
        if nargout > 1
            am = pvl_relativeairmass(90-SP.El).*Patm/101325;
            [C.CSGHI,C.CSBNI,C.CSDHI] = pvlmod_clearsky_ineichen(Loc,t,TL,SP.El,am,C.ENI);
        end
        return; 
    end
    
    lat = Loc.latitude;

    % Resample to one-minute time-steps (without having to recalculate SPA), by interpolating 
    % linearly on declination and hour-angle.
    [tt,W] = resamplestructure(datenum(t),k,'-centered');
    tt = datetime(tt,'convertfrom','datenum','TimeZone','UTC');
    p = W*double(Patm);
    ta = W*double(Ta);
    w = atan2d(W*double(sind(SP.w)),W*double(cosd(SP.w)));
    [~,el] = solarposition.eq2hor(lat,W*SP.dec,w);
    el = el + solarposition.refraction(el,p,ta);

    % Evaluate clear-sky model @ 1-min resolution
    am = pvl_relativeairmass(90-el).*p/101325;
    if ~isempty(TL), TL = W*double(TL); end
    [CSGHI,CSBNI,CSDHI] = pvlmod_clearsky_ineichen(Loc,tt,TL,el,am,W*double(C.ENI));
    C.CSGHI = single(avgdownsample(CSGHI,k));
    C.CSBNI = single(avgdownsample(CSBNI,k));
    C.CSDHI = single(avgdownsample(CSDHI,k));
    
    day = C.CSBNI > 0 & (C.CSGHI - C.CSDHI) > 0;
    
    % Correct effective-apparent-elevation to reflect CSBHI/CSBNI
    el = asind((C.CSGHI(day) - C.CSDHI(day))./C.CSBNI(day));

    % Get corresponding hour angle, and azimuth, assuming dec(el) ~ dec(a)
    a = el - solarposition.refraction(el,Patm(day),Ta(day),'-app'); % true solar position
    w = (sind(a)-sind(SP.dec(day)).*sind(lat))./(cosd(lat).*cosd(SP.dec(day))); % cos(w)
    w = max(-1,min(1,w)); % numerical precision errors can cause imaginary components otherwise
    w = sign(SP.w(day)).*acosd(w);
    az = solarposition.eq2hor(lat,SP.dec(day),w,'N2E');
 
    SP.Az(day) = az;
    SP.El(day) = el;
    SP.w(day) = w;
    
    SP.Az = single(SP.Az);
    SP.El = single(SP.El);
    SP.w = single(SP.w);
    SP.dec = single(SP.dec);
end
    
