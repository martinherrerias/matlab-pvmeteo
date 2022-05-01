function Defaults = DefaultOptions()
% Returns a structure with default simulation options
% TODO: should be replaced by a language-neutral configuration file (XML, YAML, ...)

    % Inherit defaults from matlab-utils
    PATH = './matlab-utils/general';
    Defaults = DefaultsFrom(fileparts(mfilename('fullpath')),PATH);

    Defaults.albedo = 0.2;      % Default albedo
    Defaults.soiling = 0.0;     % Default soiling factor
    Defaults.stickler = false;  % case sensitive name matching
    Defaults.useweb = true;     % Use www.geonames.org and www.soda-pro.com to retrieve information
                                % See USERACCOUNTS, GEONAMES, WPS_CAMS_MCCLEAR, and WPS_MERRA2
    
    Defaults.meteo.separationmodel = '';  % "best" given inputs (See DIFFUSE_FRACTION) 

    Defaults.minGHI = 2.0;          % Min. GHI (W/m²) for MeteoData ~obj.dark
    Defaults.minSunEl = 0.0;        % Min. solar elevation (degrees >= 0) for ~obj.dark
    Defaults.meteo.fillgaps = 0;    % Interpolate and keep(*) missing irradiance values
                                    % up to X samples ahead/before non-flagged data.
                                    % NOTE: missing ancilliary variables (e.g. Ta, Patm,..) are
                                    % interpolated in any case.
                                   
    Defaults.meteo.fitclearsky = 0;    % EXPERIMENTAL: attempt to correct bias in clear-sky model   
    Defaults.meteo.clearskycheck = 0;  % EXPERIMENTAL: flag values based on clear-sky model
                                         
    Defaults.meteo.minCSfraction = 0.05;   % Minimum required fraction of clear-sky-samples
    
    % METEOQC
    Defaults.meteoQC.basic = false;
    Defaults.meteoQC.independent = false;
    Defaults.meteoQC.sequence = {'NA','num','sunpos','weather','BSRN','CIE'};
    Defaults.meteoQC.valid = @(x) isfinite(x);
    Defaults.meteoQC.flag2nan = {};
    Defaults.meteoQC.iterations = 1;
    
    Defaults.IAM = 'martin-ruiz';  % see CHECKIAM
end

