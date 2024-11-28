classdef solarposition
% Functions for approximate solar position, mostly from Duffie & Beckman, who took many of them 
%   from Spencer (1971). Other sources quoted on the individual functions.
%
%   [1] Duffie, John A., and William A. Beckman. 2013. Solar Engineering of Thermal Processes. 
%   4th ed. Hoboken: John Wiley.
%   [2] J. W. Spencer. 1971. “Fourier Series Representation of the Position of the Sun.” 
%   Search 2 (5), 172. https://www.mail-archive.com/sundial@uni-koeln.de/msg01050.html.
%   [3] Reda, I. and Andreas, A. Solar Position Algorithm for Solar Radiaiton Applications. 
%   NREL Report No. TP-560-34302. Revised January 2008.

methods (Static)
    
[SP,C] = effectivesunpos(Loc,t,varargin)

function x = doy(t,varargin)
% X = DOY(YMD) - return Day-of-Year numbers for [year,month,day] array YMD.
% X = DOY(T) - return fractional Day-of-Year from time T, with DOY(Jan 1. 12:00 UTC) = 1.0.
% X = DOY(T,..) pass additional arguments to PARSETIME
%
%   Note that N = round(X) matches the convention of Duffie & Beckmanm (PVL_DATE2DOY)
%
% X = DOY(N) - where N is already a valid DOY [0.5 to 365.5] just returns N
%
% See also: PARSETIME, PVL_DATE2DOY

    if nargin == 1 && isnumeric(t) && ismatrix(t) && size(t,2) == 3 && ~any(mod(t,1) > 0,'all')
        x = datenum(t) - datenum([t(:,1),ones(size(t,1),2)]) + 1;
    elseif isnumeric(t) && ~any(t < 0.5 | t > 366.5,'all')
        x = t; 
    else
        t = parsetime(t,'output','datetime',varargin{:});
        t.TimeZone = 'UTC';
        x = days(t - dateshift(t,'start','year'))+0.5;
    end
end

function ds = arcdist(Ael,Aaz,Bel,Baz,r)
% ds = ARCDIST(Ael,Aaz,Bel,Baz,R)
% Returns the great-circle-arc distance between sets of points A, B located on a sphere of
% radius R. The coordinates of each point are provided in DEGREES as elevation and azimuth 
% (Xel,Xaz). If radius is omitted, default is set to 1. Examples:
%
%   ARCDIST(0,45,90,0)*180/pi - angle between a point in the plane and the zenith (90°)
%   ARCDIST(48.78,9.17,19.43,-99.13,6370) - approx. distance (km) from Stuttgart to Mexico

    if nargin < 5, r = 1; end
    compatiblesize(Ael,Aaz,Bel,Baz,r);
    Ael = double(Ael); cA = cosd(Ael); sA = sind(Ael);
    Bel = double(Bel); cB = cosd(Bel); sB = sind(Bel);
    D = double(Aaz) - double(Baz);
    cD = cosd(D); sD = sind(D);
    r = double(r);

    % ds = r.*acos(sA.*sB + cA.*cB.*cD);

    % Vicenty formula (numerically more stable)
    ds = r.*atan2(hypot(cB.*sD,cA.*sB - sA.*cB.*cD),sA.*sB + cA.*cB.*cD);
end

function CIA = cos_incidence_angle(sunel, sunaz, surftilt, surfaz)
% CIA = COS_INCIDENCE_ANGLE(SUNEL,SUNAZ,SURFTILT,SURFAZ) - all angles in degrees,
%       any consistent azimuth convention.

    CIA = sind(sunel).*cosd(surftilt)+cosd(sunel).*sind(surftilt).*cosd(sunaz-surfaz);
end

function decl = declination(doy)
% DECL = DECLINATION(DOY) - return declination (degrees) from fractional day of year
%   DOY(Jan 1. 12:00) = 1.0. Using custom fit to SPA algorithm (1990 to 2050).
% See also: DOY, HOURANGLE, ZENITH

    assert(isnumeric(doy) && ~any(doy > 366.5 | doy < 0.5,'all'),'Ivalid DOY');

        % % From Spencer (1971), in Duffie-Beckman 3rd Ed.
        % b = floor(doy)*2*pi/365; 
        % decl = (0.006918-0.399912*cos(b)+0.070257*sin(b)-0.006758*cos(2*b)+...
        %     0.000907*sin(2*b)-0.002679*cos(3*b)+0.00148*sin(3*b))*180/pi;

    % Fit to SPA (1990 to 2050), 0.09° RMS vs 0.15° of Spencer, 0.11° for pvl_ephemeris
    b = doy*2*pi/365.25;
    decl = (0.3745169-22.9084*cos(b)+4.038084*sin(b)-0.3825483*cos(2*b)+...
            0.04167191*sin(2*b)-0.1512856*cos(3*b)+0.08449235*sin(3*b));
end

function [az,el] = eq2hor(lat,decl,w,azconv)
% [AZ,EL] = EQ2HOR(LAT,DECL,W,AZCONV) - Equatorial to local-horizontal coordinate change.
%
% INPUT:
%   LAT: latitude, DECL: declination, W: hour angle. All angles in decimal degrees.
%       They can be arrays of any element-wise compatible size.
%   AZCONV: azimuth convention (e.g. 'N2E','S2W','eq90'..) default is N2E.
%
% OUTPUT:
%   AZ: azimuth, in provided convention. EL: elevation angle.
%
% REF: Duffie-Beckman 4th Ed. eq. 1.6.6 (pg. 15).
%
% See also: HOR2EQ, FIXAZIMUTH

    if nargin < 4, azconv = 'N2E'; end

    compatiblesize(lat,decl,w);
    w = solarposition.fixtoplusminus180(w);

    z = acosd(cosd(lat).*cosd(decl).*cosd(w)+sind(lat).*sind(decl));
    az = (cosd(z).*sind(lat)-sind(decl))./(sind(z).*cosd(lat));
    az = sign(w).*abs(acosd(az)); % S2W convention: S = 0°, W = 90°

    f = cosd(lat) == 0 & true(size(az));
    if any(f), az(f) = w(f); end

    az(w == 0 & decl <= lat) = 0;
    az(w == 0 & decl > lat) = 180;

    az = solarposition.fixazimuth(az,'S2W',azconv,lat);
    el = solarposition.fixtoplusminus180(90-z);
end

function [w,dec] = hor2eq(lat,az,el,azconv)
% [W,D] = HOR2EQ(LAT,AZ,EL,AZCONV) - local-horizontal to equatorial coordinate change.
%
% INPUT:
%   LAT: latitude, AZ: azimuth, in specified convention. EL: elevation angle. 
%       All angles in decimal degrees, arrays of any element-wise compatible size.
%   AZCONV: azimuth convention (e.g. 'N2E','S2W','eq90'..) default is N2E.
%
% OUTPUT:
%    DECL: declination, W: hour angle.
%
% REF: F. Vincent, Positional Astronomy, Nov-2003. [Accessed: 29-Oct-2019].
% Available: http://star-www.st-and.ac.uk/~fv/webnotes/index.html. 
%
% See also: EQ2HOR, FIXAZIMUTH

    if nargin < 4, azconv = 'N2E'; end
    az = solarposition.fixazimuth(az,azconv,'N2E',lat);

    [lat,az,el] = compatiblesize(lat,az,el);

    dec = asind(sind(el).*sind(lat) + cosd(el).*cosd(lat).*cosd(az));
    w = atan2d(-sind(az).*cosd(el),(sind(el) - sind(dec).*sind(lat))./cosd(lat));

    f = cosd(lat) == 0;
    if any(f), w(f) = az(f); end
end

function z = zenith(lat,decl,w)
    [~,el] = solarposition.eq2hor(lat,decl,w);
    z = solarposition.fixtoplusminus180(90-el);
end

function az = azimuth(lat,decl,w)
    [az,~] = solarposition.eq2hor(lat,decl,w);
end

function r = refraction(a,varargin)
% R = REFRACTION(EL,PATM,TA) - Approx. refraction angle (degrees) as a function of true 
%   solar elevation EL (degrees), atmospheric pressure PATM (Pa), and ambient temperature
%   TA (°C). 
%   NOTE: R will return the refraction R0 at apparent sunrise/sunset for points EL < -R0
%
% R = REFRACTION(A,..,'-app') - Read A as apparent-elevation angles, i.e. returning R
%   so that R = REFRACTION(A,..,'-app') ~ REFRACTION(A-R,..) within ±8"
%
% REF: J. Meeus, Astronomical algorithms, 2nd ed. Willmann-Bell, 1998. pp.106-107.

    [opt,varargin] = getflagoptions(varargin,{'-app'});
    varargin(end+1:2) = {[]};
    [Patm,Ta] = deal(varargin{:});
    if isempty(Patm), Patm = 1.01e5; end
    if isempty(Ta), Ta = 10; end
    assert(isnumeric(Patm) && isnumeric(Ta),'Unrecognized arguments');
    [a,Patm,Ta] = compatiblesize(a,Patm,Ta);

    if opt.app
    % Bennett's formula for apparent elevation 
        r90 = 1./tand(90 + (7.31./(94.4))); 
        r90 = r90-0.06*sind(14.7*r90+13);    % correction so that r(90°) = 0
        a = max(a,0);
        r = 1./tand(a + (7.31./(a+4.4)));
        r = r-0.06*sind(14.7*r+13)-r90*a/90;
    else
    % Sæmundsson's formula for true elevation
        r90 = 1.02./tand(90 + (10.3./(95.11))); % correction so that r(90°) = 0
        a0 = -0.57391415785;   % the point x for which r(-x) = x, i.e. x+r(x) = 0
        a = max(a,a0);
        r = 1.02./tand(a + (10.3./(a+5.11)))-r90*a/90; % arc minutes 
    end
    r = Patm/1.01e5.*283./(273 + Ta).*r/60;

    % 108.38/101.*238./(273-40)*0.57391415785
end

function fixedangle = fixtoplusminus180(angle)
% F = FIXTOPLUSMINUS180(A)
% Limit any angle(s) to values between -180° and 180°
    fixedangle = rem(rem(angle,360) + 540,360)-180;
end

function w0 = sunsetangle(lat, decl, varargin)
% W0 = SUNSETANGLE(LAT,DECL[,TILT]) - sunset hour angle, Duffie-Beckman 3rd Ed. pg. 16-17.
% W0 = SUNSETANGLE(LAT,DECL,'-app',PATM,TA) - refraction-corrected sunset angle (TILT = 0)
% W0 = SUNSETANGLE(..,'-finite') - Unless this flag is set, W0 will be NaN if there is no 
%   actual sunset/sunrise (arctic regions). With '-finite', dw will be 180° during arctic
%   summer (no sunset), and 0° on arctic winter (no sunrise).

    [opt,varargin] = getflagoptions(varargin,{'-finite','-app'});
    switch numel(varargin)
    case 0, tilt = 0;
    case 1, tilt = varargin{1};
    case 2
        tilt = 0;
        [Patm,Ta] = deal(varargin{2:3});
    otherwise
        error('Unrecognized syntax')
    end

    compatiblesize(lat,decl,tilt);
    arctic = any(abs(lat) > 65); % consider the possibility of refraction                

    w0 = -tand(lat).*tand(decl);
    if ~arctic
        w0 = acosd(w0);
    else
        nosunrise = w0 > 1;
        nosunset = w0 < -1;
        regular = ~(nosunrise | nosunset);
        w0(regular) = acosd(w0(regular));
        if opt.finite
            nosunriseval = 0; nosunsetval = 180;
        else
            nosunriseval = NaN; nosunsetval = NaN;
        end
        w0(nosunrise) = nosunriseval; 
        w0(nosunset) = nosunsetval; 
    end

    if ~opt.app 
        if any(tilt ~= 0)
            % get sunset angle with an artificial latitude of LAT − B"
            wt = solarposition.sunsetangle(lat - sign(lat).*tilt, decl,0);
            w0 = min(w0, wt);
        end
        return;
    end

    A = cosd(lat).*cosd(decl);
    B = sind(lat).*sind(decl);
    if arctic
        nosunrise = appsunel(0) < 0;
        nosunset = appsunel(180) > 0;
        appregular = ~(nosunrise | nosunset);

        if ~any(appregular)
            w0(nosunrise) = nosunriseval; 
            w0(nosunset) = nosunsetval; 
            return; 
        end

        arctic = ~all(appregular);
        if arctic
           w0(appregular & ~regular) = 90; % apparent sunrise/sunset only (seed §)

           A = A(appregular); B = B(appregular); w0 = w0(appregular);
           if numel(Patm) > 1
               Patm = compatiblesize(Patm,appregular); Patm = Patm(appregular); 
           end
           if numel(Ta) > 1
               Ta = compatiblesize(Ta,appregular); Ta = Ta(appregular); 
           end
        end
    end

    n = numel(w0);
    searchopt = optimoptions('lsqnonlin','Jacobpattern',speye(n),...
        'StepTolerance',1e-4,'FunctionTolerance',1e-12,'Display','none');
    w0 = lsqnonlin(@appsunel,w0,zeros(n,1),repmat(180,n,1),searchopt); %(§)

    if arctic
        w0 = revertfilter(w0,regular,[],nosunriseval);
        w0(nosunset) = nosunsetval;
    end

    function appEl = appsunel(w)
        appEl = 90 - acosd(A.*cosd(w)+B);
        appEl = appEl + solarposition.refraction(appEl,Patm,Ta);
    end
end

 function [t,dec] = solarnoon(Loc,YMD,varargin)
% T = SOLARNOON(LOC,YMD) - Calculate solar noon time (DATETIME T) for the N dates in Nx3
%   YMD [year, month, day]. Uses an internal call to PVLMOD_EPHEMERIS.
%
% T = SOLARNOON(LOC,YMD,'approx') - return approximation (±1 min), similar to Spencer 1971
% T = SOLARNOON(LOC,YMD,[delta_t],'spa') - Use PVLMOD_SPA
%
%  TODO: just replace by part of the EPHEMERIS/SPA algorithm.

    if nargin < 3 || isempty(varargin) || ~ischar(varargin{end}), mdl = 'ephemeris';
    else
        mdl = varargin{end};
        varargin(end) = [];
    end

    assert(isnumeric(YMD) && size(YMD,2) == 3,'Expecting N×3 matrix');
    % t = datetime(YMD,'TimeZone','UTC') + 0.5 - round(Loc.longitude/15)/24; % ~ noon

    doy = solarposition.doy(YMD);
    b = (doy-0.5)*2*pi/365;

    % Spencer (1971) approximation to within ±1 min
    % E = 229.2*(0.000075+0.001868*cos(b)-0.032077*sin(b)-0.014615*cos(2*b)-0.04089*sin(2*b));

    % Custom fit to SPA
    E = 0.011990+0.534137*cos(b)-7.325240*sin(b)-3.484400*cos(2*b)-9.292518*sin(2*b)...
        -0.091508*cos(3*b)-0.328605*sin(3*b)-0.141505*cos(4*b)-0.179312*sin(4*b);

    t = datetime(YMD,'TimeZone','UTC') + hours(12) - minutes(4*Loc.longitude + E); % ~ noon

    switch lower(mdl)
    case 'approx'
        if nargout > 1
            doy = mod(doy - (4*Loc.longitude + E)/1440 - 0.5,366)+0.5; 
            dec = solarposition.declination(doy);
        end
        return;
    case 'ephemeris'
        [~,~,~,w,dec] = pvlmod_ephemeris(t,Loc);
    case 'spa'
        % Use hour angle from SPA to reduce w to SPA ± 6e-5 (SPA precision is ±3e-4)
        [~,~,~,w,dec] = pvlmod_spa(t,Loc,varargin{:});
    otherwise
        error('Unknown option %s',mdl);
    end
    t = t - hours(w/15);
 end

function [dw,t,sr,ss] = daytime(Loc,YMD,varargin)
% [H,SN,SR,SS] = DAYTIME(LOC,YMD)
% [..] = DAYTIME(LOC,YMD,'-app',[PRESSURE,TEMPERATURE])
%   Calculate [apparent] daylight hours (H), solar noon (SN), sunrise (SR) and
%   sunset (SS) times for the N dates in N×3 matrix YMD
% NOTE: use SOLARPOSITION.SUNSETANGLE if all you need is approx. H 
% See also: SOLARPOSITION.SOLARNOON, SOLARPOSITION.SUNSETANGLE

    % remove flags, keep pressure and temperature, if available
    [~,optargs] = getflagoptions(varargin,{'-finite','-app'});

    [t,dec] = solarposition.solarnoon(Loc,YMD,optargs{:});
    dw = solarposition.sunsetangle(Loc.latitude,dec,varargin{:});
    dw = dw/15;
    if nargout > 2, sr = t - hours(dw); end
    if nargout > 3, ss = t + hours(dw); end
end
function sr = sunrisetime(varargin)
% SR = SUNRISETIME(LOC,YMD,..) - alias for SOLARPOSITION.DAYTIME
    [~,~,sr] = daytime(varargin{:}); end
function ss = sunsettime(varargin)
% SS = SUNSETTIME(LOC,YMD,..) - alias for SOLARPOSITION.DAYTIME
    [~,~,~,ss] = daytime(varargin{:}); end

function w = hourangle(t,longitude, varargin)
% W = HOURANGLE(DOY,LON,TIMEZONE,TIME)
% From Spencer (1971)-Duffie-Beckman 4th Ed. 1.5.
% DOY: fractional day of year (1 to 366)
% LON: longitue (degrees)
% TIMEZONE: time zone (hours)
% TIME: optional, fractional days, otherwise calculated from fractional days of DOY

    t = parsetime(t,varargin{:});
    t.TimeZone = 'UTC';

    b = days(t-dateshift(t,'start','year'))*2*pi/365;
    E = 229.2*(0.000075+0.001868*cos(b)-0.032077*sin(b)-0.014615*cos(2*b)-0.04089*sin(2*b));

    solartime = hours(t-dateshift(t,'start','day'))-(longitude*4+E)/60;
    w = (solartime-12)*15;
end

function newaz = fixazimuth(az,orgconv,endconv,lat)
% newaz = FIXAZIMUTH(az,orgconv,endconv,lat)
%   Turn azimuth angle(s) az in convention 'orgconv' into convention 'endconv'.
%
%   az: decimal degree angle(s).
%   orgconv, endconv: azimuth convention 3-char strings of the form 'A2B', where A and B  
%       are cardinal-points {N,E,S,W}, eg. N2E means N = 0°, E = 90°, W = 270° i.e. the
%       astronomical convention, whereas S2W is also common in the solar-energy sector.
%
%       Strings of the form 'eqX' or simply a numeric X are read as decimal degree angles CCW 
%       from equator, e.g. eq0 / eq90 are equal to S2E / W2S if lat >= 0 or N2W / E2N if lat < 0.
%       A fixed-mount solar array will typically have az ~ 0 in eq0 convention.
%
%   lat: latitude, is required if either convention is numeric (offset from equator)
%
% [newaz,endconv] = FIXAZIMUTH(az,orgconv,endconv,lat,[default]) - empty orgconv/endconv are
%   assumed to be = default (SimOptions.AzConv, if missing), with a warning.
%
% See also: solarposition.fixtoplusminus180

    if nargin < 4, lat = []; end
    
    [refangle(1),dirsign(1)] = sortoutconvention(orgconv,lat);
    [refangle(2),dirsign(2)] = sortoutconvention(endconv,lat);

    newaz = prod(dirsign)*az + dirsign(2)*diff(refangle);
    newaz = solarposition.fixtoplusminus180(newaz);

    function [refangle,direction] = sortoutconvention(conv,lat)
    % Find meaning from a convention string as described above.
    % refangle: azimuth in N2E convention for zero-azimuth cardinal point
    % direction: 1 for clockwise (NWSE..), -1 for counter-clockwise (NESW..)

        if isnumeric(conv)
            validateattributes(conv,{'numeric'},{'scalar','real','<=',360,'>=',-360},'','conv');
            assert(~isempty(lat),'fixazimuth:eq','Latitude required for eq+ convention');
            validateattributes(lat,{'numeric'},{'scalar','real','<=',90,'>=',-90},'','lat');
            refangle = conv + 180*(lat >= 0); direction = 1;
        else
            validateattributes(conv,{'char'},{'nonempty','size',[1,NaN]'},'','conv');
            switch upper(conv)
            case 'N2E', refangle = 0; direction = -1;
            case 'E2S', refangle = 90; direction = -1;
            case 'S2W', refangle = 180; direction = -1;
            case 'W2N', refangle = 270; direction = -1;
            case 'N2W', refangle = 0; direction = 1;
            case 'W2S', refangle = 270; direction = 1;
            case 'S2E', refangle = 180; direction = 1;
            case 'E2N', refangle = 90; direction = 1;
            otherwise
                x = textscan(upper(conv),'EQ%f');
                if ~isempty(x{1}) && x{1} <= 360 && x{1} >= -360
                    [refangle,direction] = sortoutconvention(x{1},lat);
                else
                    error('fixazimuth:A2B',['Invalid azimuth convention string. Expecting ',...
                          '''A2B'', where A, B are cardinal points {N,S,E,W} 90° appart']);
                end
            end
        end
    end
end

function C = cardinal(az,tol)
% C = CARDINAL(AZ,[TOL]) - Return a cell array of (tertiary inter-) cardinal directions for
%   azimuth angles AZ (N2E convention).
%
% EXAMPLE:
%   az = 0:360/16:359;
%   plothalfdome();
%   text(sind(az)*1.1,cosd(az)*1.1,solarposition.cardinal(az));
    
    if nargin < 2, tol = 360/64; end
    validateattributes(tol,'numeric',{'scalar','real','positive','<=',360/64});
    
   L = {'S' 'SbW' 'SSW' 'SWbS' 'SW' 'SWbW' 'WSW' 'WbS' ...
        'W' 'WbN' 'WNW' 'NWbW' 'NW' 'NWbN' 'NNW' 'NbW' ...
        'N' 'NbE' 'NNE' 'NEbN' 'NE' 'NEbE' 'ENE' 'EbN' ...
        'E' 'EbS' 'ESE' 'SEbE' 'SE' 'SEbS' 'SSE' 'SbE'};
   az = solarposition.fixtoplusminus180(az);
   [match,idx] = ismembertol(az,-180:360/32:180-tol/2,tol,'datascale',1);
   
   C = cell(size(az));
   C(match) = L(idx(match));
   C(~match) = arrayfun(@(x) num2str(x,'%0.0f°'),az(~match),'unif',0);
end

function [az,el,w,d,T] = sunposgrid(lat,varargin)
% [AZ,EL,W,D,T] = SUNPOSGRID(LAT,LATTICE) - Generates a triangular grid of solar-positions 
%   for latitude LAT, evenly spaced by a given LATTICE.
%
%   Points are generated over a rectangular grid of hour-angles and declinations, and
%   rows (constant declination) are alternatively offset by ±LATTICE/4 to create a 
%   triangular mesh. Points on the grid that will never be used for interpolation above
%   the horizon (i.e. points not connected on a triangulation* to any other points above
%   the horizon) are finally removed.
%
% INPUT: LAT (degrees north), LATTICE (positive scalar, degrees), or 2-vector specifying
%   grid spacing on hour-angle and declination, respectively.
% OUTPUT: azimuth AZ (N2E), elevation EL, hour-Angle W, and declination D. All  angles in 
%   degrees. T triangulation used for filtering.
%
% SUNPOSGRID(..,'minSunEl',E0) - Use an alternative limit (degrees above/below horizon)
%   to filter out useful points.
%
% SUNPOSGRID(..,'-rectangular') - Don't apply ±LATTICE/4 offset to rows. 
% SUNPOSGRID(..,'-matrix') - In combination with '-rectangular' flag, returns AZ, EL
%   angles as matrices (rows of equal declination, columns of equal hour-angle). Useless
%   values are flagged as NaN's instead of removed. W, D are reduced to a row- and column-
%   vector, respectively.
% 
% SUNPOSGRID(..,'-plot') - Generate a plot of the grid on GUIfigure('sunposgrid')
%
% See also: SUNPOSMESH

    % declination steps shouldn't go beyond ±23.5°, but a small margin allows meteo-data 
    % to have a slightly different latitude
    DEF.lattice = 5;
    DEF.azconv = 'N2E';
    DEF.maxdec = 24.0;
    DEF.minSunEl = 0;

    DEF.rectangular = false;
    DEF.matrix = false;
    DEF.plot = false;

    opt = parseoptions(varargin,{'-rectangular','-matrix','-plot'},DEF,'dealrest',1);
    % [opt,varargin] = getflagoptions(varargin,{'-rectangular','-matrix','-plot'});
    % opt = getpairedoptions(varargin,completeoptions(opt,DEF),'dealrest',1);

    parsestruct(opt,{'lattice'},'-n','-r','-f','-v',@(x) any(numel(x) == [1,2]));
    parsestruct(opt,{'maxdec','minSunEl'},'-n','-r','-f','-s');
    parsestruct(opt,{'rectangular','matrix','plot'},'-s','-l');
    if isscalar(opt.lattice)
       if ~opt.rectangular, opt.lattice(2) = sqrt(3)/2*opt.lattice; 
       else, opt.lattice(2) = opt.lattice;
       end 
    end

    wstep = opt.lattice(1);
    decstep = opt.lattice(2);

    d = ceil(2*opt.maxdec/decstep)*decstep; 
    d = -d/2:decstep:d/2;

    % hour angles for z = minSunEl° - max(lattice) < 0
    z0 = opt.minSunEl - max(decstep,wstep);
    w0 = (sind(z0)-sind(lat)*sind(d))./(cosd(lat).*cosd(d));
    w0 = acosd(max(min(w0,1),-1));

    w = ceil(2*max(w0)/wstep)*wstep;
    % if ~opt.rectangular, w = w + wstep; end
    w = -w/2:wstep:w/2;

    [w,d] = meshgrid(w,d);
    if ~opt.rectangular
    % Offset alternate rows by ±dw/4
        w = w + wstep/4 - wstep/2*mod(1:numel(w0),2)';
    end

    % Get azimuth and elevation for existing equatorial coordinates
    [az,el] = solarposition.eq2hor(lat,d(:),w(:),opt.azconv);

    s = [cosd(az).*cosd(el),sind(az).*cosd(el),sind(el)];
    T = convhull(s);

    % Remove triangles with edges longer than reasonable
    if ~opt.rectangular, d0 = wstep; else, d0 = hypot(wstep,decstep); end
    d0 = 1.2*d0*pi/180;
    edgeok = @(j,k) rssq(s(T(:,j),:) - s(T(:,k),:),2) <= d0;
    T = T(edgeok(1,2) & edgeok(2,3) & edgeok(1,3),:);

    f = el > opt.minSunEl;  % vertices with elevation > minSunEl
    T = T(any(f(T),2),:);   % facets with at least one point > minSunEl
    [f,~,idx] = unique(T);   % vertices connected to at least one point > minSunEl

    if opt.matrix && opt.rectangular
        f = full(sparse(f,1,true,numel(az),1));
        az = reshape(az,size(w)); el = reshape(el,size(w));
        az(~f) = NaN; az(~f) = NaN;
        w = w(1,:); 
        d = d(:,1);
    else
        w = w(f); d = d(f);
        az = az(f); el = el(f);
        T = reshape(idx,[],3);
    end

    if opt.plot
        GUIfigure('sunposgrid'); clf(); hold on;
        plothalfdome();
        s = [cosd(az).*cosd(el),sind(az).*cosd(el),sind(el)];
        scatter3(s(:,1),s(:,2),s(:,3),2,'r');
        patch('Faces',T,'Vertices',s,'facecolor','none','edgecolor',[1 0.8 0.8],'linestyle',':');
    end 
end

function [az,el,w] = meanposition(lat,w,dec,Patm,Ta,dt,varargin)
% [AZ0,EL0,W0] = MEANPOSITION(LAT,W,DEC,PATM,TA,DT,AZCONV) - Calculate an approx. "mean" 
%   position of the sun within finite intervals (W(j)-dW/2,W(j)+dW/2), where dW is the
%   angular displacement of the sun within successive time-intervals.
%
%   ATENTION: while better than the center-of-interval approximation, MEANPOSITION
%   corresponds to a Top-of-Atmosphere-weighted position, which is worse than a clear-
%   sky-model-weighted position (see EFFECTIVESOLARPOSITION and Blanc & Wald, 2016).
%
%   The mean solar-elevation angle EL0(j) corresponds to the point in the solar path
%   W(j)-dW/2 -> W(j)+dW/2 for which sin(EL0(j)) is representative of the whole integrated
%   interval:
%                    1  / w(j)+dW/2
%      sin(EL0(j)) = -- |   max(0,sin(el))dw
%                    dW / w(j)-dW/2
%
%   For small time-steps (~1 min) and intervals that don't include the sunset/sunrise, the
%   approximation EL0 ~ el(j) might be sufficient, but in other cases it can produce errors
%   of order O(1e-2).
%
%   MEANPOSITION uses the analytic solution to the integral along a solar path: 
%
%       sin(el) ~ A·cos(W) + B 
%
%   Without refraction, A = cos(LAT)·cos(DEC) and B = sin(LAT)·sin(DEC) are functions of 
%   declination and latitude only; due to zenith-dependent atmospheric refraction near
%   the horizon, B must be adjusted, and an error of O(1e-5) on sin(EL0) should still be
%   expected.
%
% INPUT:
%
%   lAT - location latitude (degrees north)
%   W - N-vector of solar hour angles (degrees) at the center of each interval.
%   DEC - N-vector of declination angles (degrees) at the center of each interval
%   PATM, TA - (optional) N-vectors of atmospheric pressure (Pa) and temperature (°C),
%       to account for refraction near the horizon. Defaults are 1.01e5 Pa and 10°C.
%   DT - (optional) time-step duration (minutes), used to calculate dW. The default is to 
%       use dW = mode(diff(W)) directly. Note this doesn't make sense for non-uniform
%       time-steps.
%   AZCONV - (optional) Azimuth convention key (see SOLARPOSITION.FIXAZIMUTH)
%
% OUTPUT: [AZ0,EL0,W0] - N-vectors of effective azimuth, elevation, and hour angles 
%   (degrees). Azimuth convention is set by AZCONV or the default in FIXAZIMUTH.
%
% REFERENCES:
%   Blanc, P., Wald, L., 2016. On the effective solar zenith and azimuth angles to use 
%   with measurements of hourly irradiation. Advances in Science and Research 13, 1–6.
%   https://doi.org/10.5194/asr-13-1-2016
%
% See also: SOLARPOSITION.EQ2HOR, EFFECTIVESOLARPOSITION

    narginchk(3,7);

    if nargin < 4 || isempty(Patm), Patm = 1.01e5; end
    if nargin < 5 || isempty(Ta), Ta = 10; end

    [w,dec,Patm,Ta] = compatiblesize(w,dec,Patm,Ta);

    if nargin < 6 || isempty(dt)
        assert(numel(w) > 1,'Cannot guess dT with numel(W) < 2');
        dw = diff(w);
        dw(dw < 0 | dw > 180) = [];
        dw = mode(dw)*pi/180; % integration interval (radians)
        dt = dw*1440/(2*pi);
    else
        dw = dt*2*pi/1440;
    end
    if ~any(round(dt) == [1,2,2.5,5,10,15,20,30,45,60])
        warning('effectivesp:common','dT = %0.1f min is not a common integration interval',dt)
    elseif dt < 0.5
        warning('effectivesp:small','for dT = %0.1f s you might as well use the center of the interval',dt*60)
    end

    e0 = solarposition.refraction(-0.5739,Patm,Ta); % ~ refraction at app. horizon
    for j = 1:10
        err = solarposition.refraction(-e0,Patm,Ta)-e0;
        e0 = e0 + err;
        % GUIfigure('debug'); clf(); hist(err);
        if all(abs(err) < 1e-6), break; end
    end

    A = cosd(lat)*cosd(dec);
    B = sind(lat)*sind(dec);
    ws = acos((sind(-e0)-B)./A); % sunset/sunrise hour angle, including refraction of ~34'

    wB = max(-ws,min(ws,w*pi/180 + [-1,1]*dw/2)); % actual integration limits (radians)

    % filter out any points completely below the horizon
    day = wB(:,1) ~= wB(:,2);
    wB = wB(day,:);
    a = A(day);
    b = B(day);

    w0 = mean(wB,2);               % at the center of the integration interval...
    el = asind(a.*cos(w0) + b);    % true elevation
    el = el + solarposition.refraction(el,Patm(day),Ta(day));
    B_mod = sind(el) - a.*cos(w0); % B, adjusted to include refraction

    C = (a.*diff(sin(wB),1,2)+B_mod.*diff(wB,1,2))/dw; % mean cos(z)

    % Update w0 as the mean of w weighted by cos(z)
    w0 = (a.*diff(wB.*sin(wB)+cos(wB),1,2)+B_mod/2.*diff(wB.^2,1,2))./(C*dw);

    % NOTE: tempting to iterate (recalculate B_mod @ w0) but actually lowers precision

    % For night points,use apparent elevation at the center of interval
    el = asind(A.*cosd(w) + B)+e0;
    [az,~] = solarposition.eq2hor(lat,dec,w,varargin{:});

    % For day points, replace by effective apparent-elevation
    el(day) = asind(C);
    [az(day),~] = solarposition.eq2hor(lat,dec(day),w0*180/pi,varargin{:});

    w(day) = w0;
end
function n = mdoy(m)
% N = MDOY(M) - Recommended representative days for Month(s) M.
% N = MDOY() - equivalent to MDOY(1:12)
% from Klein (1977)-Duffie Beckman 3rd Ed.

    if nargin < 1, m = 1:12; end

    n0 = [17, 47, 75, 105, 135, 162, 198, 228, 258, 288, 318, 344];
    n = n0(m);
    if isvector(m) && size(m,1) > 1, n = n'; end
end

function e = extraradiation(t,C)
% E = EXTRARADIATION(T,C) - return extraterrestrial normal irradiance for times T and
%   solar constant C (default 1361 W/m², base on Kopp & Lean, 2011). Uses the eqn. form  
%   of Spencer (1971) but custom-fit to SPA calculated sun-earth distance factor 1/r², 
%   reducing the error by ~0.6 W/m², (which, btw. matters much less than the change of
%   solar constant from the default 1367 W/m²).
%
%  [1] Kopp, G., Lean, J.L., 2011. A new, lower value of total solar irradiance: 
%   Evidence and climate significance. Geophysical Research Letters 38.

    if nargin < 2, C = 1361.6; end % 1360.8 min. + 1.6/2 of 11-year cycle variation [1]

    % Spencer (1971) and Duffie-Beckman 3rd Ed. 1.4 (pg. 9)
    % b = floor(n)*2*pi/365;
    % e = 1367*(1.00011+3.4221e-2*cos(b)+1.28e-3*sin(b)+7.19e-4*cos(2*b)+7.7e-5*sin(2*b));

    b = solarposition.doy(t)*2*pi/365;
    e = C*(1.00014 + 3.33513e-02.*cos(b)+1.92329e-03.*sin(b)+...
           6.93126e-04.*cos(2*b)+7.81281e-05.*sin(2*b));
end

end
end

% function HayDavies(DHI, DNI, AOI, zenith, tilt, Optional Ea = 1367)
% 
%     Dim Ai  %Anisotropy index
%     Dim Rb  %Ratio of beam irradiance on the tilted surface vs. on a horizontal surface
% 
%     Ai = DNI/Ea
%     Rb = max(cos(AOI*pi/180), 0)/max(cos(zenith*pi/180), 0.01745)
%     HayDavies = DHI*(Ai*Rb+(1-Ai)*(1+cos(pi*tilt/180))/2)
% 
% end
%
% function HDKR(DHI, GHI, AOI, zenith, tilt, Optional rho = 0.2, Optional Eo = 1367)
% %Gpoa = HDKR(DHI,GHI,AOI,zenith,tilt,rho,Eo)
% %Hay-Davies-Klucher-Reindl Model (Duffie-Beckman 3rd Ed. pg. 93)
% 
%     Dim Ai  %Anisotropy index
%     Dim Rb  %Ratio of beam irradiance on the tilted surface vs. on a horizontal surface
%     Dim f   %horizon brightening modulating factor
%     Dim Eh
%     
%         if GHI <= 0 
%             HDKR = 0
%             return
%         end
%         
%         Eh = Eo*cosd(zenith)
%         
%         if Eh < GHI 
%             GHI = max(GHI, DHI)
%             DHI = GHI
%             Eh = 1
%         end
%         
%         Ai = (GHI-DHI)/Eh
%         f = Sqr((GHI-DHI)/GHI)
%         Rb = max(cosd(AOI), 0)/max(cosd(zenith), 0.01745)
%         
%     HDKR = (GHI+DHI*(Ai-1))*Rb+DHI*(1-Ai)*(1+cosd(tilt))/2*(1+f*(sind(tilt/2)) ^ 3)+rho*GHI*(1-cosd(tilt))/2
% 
% end
% 
% function GroundReflected(GHI, albedo, tilt)
%     GroundReflected = GHI*albedo*(1-cos(pi*tilt/180))/2
% end
% 
% function std_atm(elev)
%     std_atm = (101325*(1-0.0000225577*elev) ^ 5.25588)*0.01
% end
% 
% function approx_monthly_diffuse_fraction(Kt, w_sunset)
% %Erbs et al. (1982), from Duffie-Beckman 4th Ed. eq. 2.12.1 (pg. 80).
% 
%     Dim DF
%     
%     if Kt > 0.8 Or Kt < 0.3 
%         approx_monthly_diffuse_fraction = Error(xlErrValue)
%     end
%         
%     if w_sunset <= 81.4 
%         DF = 1.391-3.56*Kt+4.189*Kt ^ 2-2.137*Kt ^ 3
%     else
%         DF = 1.311-3.022*Kt+3.427*Kt ^ 2-1.821*Kt ^ 3
%     end
%     
%     approx_monthly_diffuse_fraction = DF
%     
% end
% 
% function approx_hourly_diffuse_fraction(Kt)
% %Erbs et al. 1982, de Duffie-Beckman 3rd. Ed. pg. 76
% 
%     Dim Id_I
%     
%     if Kt <= 0.22 
%         Id_I = 1-0.09*Kt
%     else
%         if Kt <= 0.8 
%             Id_I = 0.9511-0.1604*Kt+4.388*Kt ^ 2-16.638*Kt ^ 3+12.336*Kt ^ 4
%         else
%             Id_I = 0.165
%         end
%     end
%     
%     approx_hourly_diffuse_fraction = Id_I
%     
% end
% 
% function Rb(lat, tilt, Ndoy)
% %Azimuth 0º, Duffie-Beckman 3rd Ed. 1.8 (pg. 24)
%     Dim b, decl, w_sunset
%     
%     b = (Ndoy-1)*2*pi/365
%     decl = approx_declination(Ndoy)*pi/180
%     w_sunset = approx_sunset_w(lat, decl, tilt)*pi/180
%     lat = lat*pi/180
%     tilt = tilt*pi/180
%     
%     Rb = (cos(lat-tilt)*cos(decl)*sin(w_sunset)+w_sunset*sin(lat-tilt)*sin(decl))/_
%          (cos(lat)*cos(decl)*sin(w_sunset)+w_sunset*sin(lat)*sin(decl))
%     
% end
% 
% function Rb_hour(lat, decl, w, w_sunset, Optional beta = 0, Optional surf_az = 0)
% %Rb = Rb_hour(latitude, declination, hour_angle, sunset_hour_angle, surf_tilt, surf_azimuth)
% %Average Rb (Incident/Horizontal Beam irradiance) over w ± 7.5° assuming constant sky transmissivity
% %Duffie-Beckman 4th Eq. 2.14.6 (pg. 88)
% 
%     Dim w1, w2
%     Dim A1, a2, a3, b
%     
%     w1 = w-15/2
%     if w1 < -w_sunset  w1 = -w_sunset
%     
%     w2 = w+15/2
%     if w2 > w_sunset  w2 = w_sunset
%     
% 
%     A1 = (sind(decl)*sind(lat)*cosd(beta)-sind(decl)*cosd(lat)*sind(beta)*cosd(surf_az))*pi/180*(w2-w1)
%     a2 = (cosd(decl)*cosd(lat)*cosd(beta)+cosd(decl)*sind(lat)*sind(beta)*cosd(surf_az))*(sind(w2)-sind(w1))
%     a3 = -(cosd(decl)*sind(beta)*sind(surf_az))*(cosd(w2)-cosd(w1))
% 
%     b = (cosd(lat)*cosd(decl))*(sind(w2)-sind(w1))+(sind(lat)*sind(decl))*pi/180*(w2-w1)
% 
%     Rb_hour = (A1+a2+a3)/b
%     
% end
% 
% function RbH_hour(lat, decl, w, w_sunset)
% %Rb = RbH_hour(latitude, declination, hour_angle, sunset_hour_angle, surf_tilt, surf_azimuth)
% %Average BNI/BHI (Normal/Horizontal Beam irradiance) over w ± 7.5° assuming constant sky transmissivity
% %Modified from Duffie-Beckman 4th Eq. 2.14.6 (pg. 88)
% 
%     Dim w1, w2
%     Dim a, b
%     
%     if abs(w) > w_sunset+7.5 
%         RbH_hour = -1
%         return
%     end
%     
%     w1 = w-15/2
%     if w1 < -w_sunset  w1 = -w_sunset
%     
%     w2 = w+15/2
%     if w2 > w_sunset  w2 = w_sunset
%     
%     a = pi/180*(w2-w1)
%     b = (cosd(lat)*cosd(decl))*(sind(w2)-sind(w1))+(sind(lat)*sind(decl))*pi/180*(w2-w1)
%     RbH_hour = a/b
%     
% end
% 
% function approx_hourly_temp(h, Tavg, Kt)
% %From: Erbs et al. 1983, "Estimation of degree-days and ambient temperature bin data from monthly average temperatures"
% %Taken from: Energy Efficient Buildings with Solar and Geothermal Resources
% 
% %h: hour, Excel std. [days]
% %Tavg: monthly average temperature (ºC)
% %Kt: monthly average clearness index
% 
%     Dim t
%     Dim a
%     Dim TH
%     
%     a = 25.8*Kt-5.23        %Good aproximation only for continental regions and mid-latitudes!
%     t = 2*pi/24*(h*24-1)
%     TH = 0.4632*cos(t-3.805)+0.0984*cos(2*t-0.36)+0.0168*cos(3*t-0.822)+0.0138*cos(4*t-3.1513)
% 
%     approx_hourly_temp = Tavg+a*TH
% end
% 
% function approx_clearsky_BNI(Alt, Rb, Optional Ea = 1367, Optional ClimType As Variant = 0)
% %Hottel (1976)formula for estimating the beam radiation transmitted through clear atmospheres
% %From Duffie-Beckman 4th ed. Eq. 2.8.1
% 
% %Alt is altitude above sea level (meters)
% %Rb is the Beam-Normal-to-Horizontal irradiance ratio, can be approximated by 1/cos(zenith_angle)
% %Ea is extraterrestrial normal irradiance (see approx. extra irradiance)
% %ClimType is a modifying variable (integer or string):
% 
% %ClimType           Climate Type         r0   r1   rk
% %0/"None"          -                  1.00 1.00 1.00
% %1/"Tropical"       Tropical            0.95 0.98 1.02
% %2/"MidLatSummer"   Midlatitude summer  0.97 0.99 1.02
% %3/"MidLatWinter"   Midlatitude winter  1.03 1.01 1.00
% %4/"SubarcSummer"   Subarctic Summer    0.99 0.99 1.01
% 
%     Dim a0, A1, k
%     Dim r0, r1, rk As Variant
%     Dim typ
%     Dim basura As String
%     
%     basura = VarType(ClimType)
%     r0 = Array(1#, 0.95, 0.97, 1.03, 0.99)
%     r1 = Array(1#, 0.98, 0.99, 1.01, 0.99)
%     rk = Array(1#, 1.02, 1.02, 1#, 1.01)
%     
%     if VarType(ClimType) = vbString 
%         Select Case ClimType
%         Case "None"
%             typ = 0
%         Case "Tropical"
%             typ = 1
%         Case "MidLatSummer"
%             typ = 2
%         Case "MidLatWinter"
%             typ = 3
%         Case "SubarcSummer"
%             typ = 4
%         Case else
%             GoTo Err
%         end Select
%     Elseif IsNumeric(ClimType) 
%         typ = ClimType
%         if Fix(typ) <> typ And (typ > 4 Or typ < 0)  GoTo Err
%     else
%         GoTo Err
%     end
%     
%     a0 = r0(typ+1)*(0.4237-0.00821*(6-Alt/1000) ^ 2)
%     A1 = r1(typ+1)*(0.5055+0.00595*(6.5-Alt/1000) ^ 2)
%     k = rk(typ+1)*(0.2711+0.01858*(2.5-Alt/1000) ^ 2)
% 
%     approx_clearsky_BNI = Ea*(a0+A1*Exp(-k*Rb))
%     return
% Err:
%     approx_clearsky_BNI = CVErr(xlErrValue)
% end
% 




