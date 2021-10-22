function [s,u] = whatis(x,varargin)
% METEODATA.WHATIS(NAME) display description for known variable NAME
% [S,U] = METEODATA.WHATIS(NAME) - return description string S and units U for NAME.
% METEODATA.WHATIS(PATTERN,'-regexpx') display all regexp pattern matches.
%
% See also: METEODATA.FINDFIELDS

    if ~ischar(x)
        x = cellstr(x);
        [s,u] = cellfun(@MeteoData.whatis,x,'unif',0);
        return;
    end

    switch x
    case 'DHI', u = 'W/m²'; s = 'diffuse horizontal irradiance';
    case 'BNI', u = 'W/m²'; s = 'direct (beam) normal irradiance';
    case 'GHI', u = 'W/m²'; s = 'global horizontal irradiance';
    case 'GTI', u = 'W/m²'; s = 'global tilted irradiance';
    case 'USW', u = 'W/m²'; s = 'upwelling short-wave irradiance (albedo)';

    case 'ENI', u = 'W/m²'; s = 'Direct normal irradiance at the top of the atmosphere (extraterrestrial).';
    case 'sunel', u = 'deg'; s = 'Refraction-corrected solar elevation angle.';
    case 'sunaz', u = 'deg'; s = 'Solar azimuth angle (N2E convention)';
    case 'hourangle', u = 'deg'; s = 'Solar hour angle';
    case 'declination', u = 'deg'; s = 'Solar declination';

    case 'Patm', u = 'Pa'; s = 'Atmospheric pressure at ground level';
    case 'RH', u = ''; s = 'Relative humidity at ground level';
    case 'tpw', u = 'mm'; s = 'atmospheric water contents';
    case 'windir', u = 'deg'; s = 'Direction of origin of wind (N2E convention).';
    case 'vw', u = 'm/s'; s = 'Wind speed measured at standard 10 meter height.';
    case 'Ta', u = '°C'; s = 'ambient dry bulb temperature.';

    case 'soiling', u = ''; s = 'irradiance loss fraction due to soiling';
    case 'CSGHI', u = 'W/m²'; s = 'Clear-Sky global horizontal irradiance';
    case 'CSBNI', u = 'W/m²'; s = 'Clear-Sky direct normal irradiance';
    case 'CSDHI', u = 'W/m²'; s = 'Clear-Sky diffuse horizontal irradiance';
    case 'clearsky', u = '0,1'; s = 'Clear-sky sample';

    case 'TL', u = ''; s = 'linke turbidity';
    case 'AMa', u = ''; s = 'absolute airmass';
    case 'AOD', u = ''; s = 'aerosol optical depth';

    case 'kd', u = ''; s = 'diffuse fraction: DHI/GHI';
    case 'kt', u = ''; s = 'clearness index: GHI / ( ENI cos(z) )';
    case 'kn', u = ''; s = 'direct clearness index: BNI / ENI';
    case 'albedo', u = ''; s = 'upwelling / GHI';

    case 'Ts', u = '°C'; s = '(generic) sensor temperature';
    case 'Tm', u = '°C'; s = 'Module temperature';

    % case '', s = 'F1';
    % case '', s = 'F2';
    % case '', s = '... other diffuse components? ';
    otherwise
        u = '?';
        s = 'Unknown variable';
        msg = [s ' "' x '"'];
        y = MeteoData.findfields(x,'-soft',varargin{:});
        if ~strcmp(x,y)
            y = cellfun(@(x) [x ', ' MeteoData.whatis(x)],cellstr(y),'unif',0);
            msg = sprintf('%s, try:\n\t%s',msg,strjoin(y,'\n\t'));
        end
        if nargout > 0
            warning('meteodata:whatis',msg);
        else
            s = msg;
        end
    end
    if ~isempty(u) && ~strcmp(u,'?')
        s = [s ' [' u ']'];
    end
end