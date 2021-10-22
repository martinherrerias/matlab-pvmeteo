function varargout = findfields(S,varargin) 
% R = METEODATA.FINDFIELDS(S) - Rename fields of S according to list of common aliases
% [R,I] = METEODATA.FINDFIELDS(S) - Return also an index I of the actual substitutions
% [M,I,U] = METEODATA.FINDFIELDS(S) - Return a partial copy M and index I with fields that 
%   have been matched/renamed, and a residual U with all other fields. 
%
% [..] = METEODATA.FINDFIELDS(S,'name','source',..,['-stickler']) - provide explicit pairs of
%   standardized* names and sources to map them from. Unless '-stickler' flag is set, names
%   are permisive, i.e. they are themselves parsed through the predefined list of common
%   aliases: stdname <- name <- alias. Otherwise case-sensitive matches are expected.
%
% [..] = METEODATA.FINDFIELDS(.. ,varargin) - provide additional arguments to RENAMEFIELDS,
%   e.g. MeteoData.findfields('clearness.index','-regexpx','-soft') will return {kt,kn}
%   whereas MeteoData.findfields('clearness index') will just crash.
%
% NOTES:
%   This is just a call to RENAMEFIELDS searching for common aliases for MeteoData properties,
%   meant to simplify compatibility with other libraries and old code versions. The list of
%   aliases should probably expand... maybe use regexp expressions by default.
%
% See also: METEODATA.WHATIS, RENAMEFIELDS, PARSELIST, METEODATA.VARNAMES

    validateattributes(S,{'MeteoData','struct','table','timetable','char','cell','string'},{});

    [opt,varargin] = getflagoptions(varargin,'-stickler');
    if opt.stickler
    % Case sensitive matches only
    
        alias = MeteoData.varnames('all');
        [alias,varargin] = getpairedoptions(varargin,alias,alias);
        if isstruct(S) || isobject(S)
            [varargout{1:2}] = renamefields(S,alias,'cat',2,varargin{:});
        else
            [varargout{1:2}] = parselist(S,alias,'-matchcase',varargin{:});
        end
        
    else
    % Everything goes... reality is fluid, man.

        ALIAS.location = {'location','loc'};
        ALIAS.info = {'info'};
        ALIAS.options = {'options','opt'};
        ALIAS.sensors = {'sensors','sens'};
        ALIAS.flags = {'flags'};

        ALIAS.t = {'t','time','date_time','timestamp','date'};
        ALIAS.timestep = {'timestep','dt'};
        ALIAS.interval = {'interval','summarization'};

        ALIAS.DHI = {'dhi'};
        ALIAS.BNI = {'bni','dni'};
        ALIAS.GHI = {'ghi'};
        ALIAS.GTI = {'gti'};

        ALIAS.ENI = {'eni','dni_toa','Ea','etrn'};
        ALIAS.sunel = {'sunel','apparent_solar_elevation'};
        ALIAS.sunaz = {'sunaz','azm','zenith'};
        ALIAS.hourangle = {'hourangle','w'};
        ALIAS.declination = {'declination','dec'};

        ALIAS.kt = {'kt','clearness_ghi','clearness_index'};
        ALIAS.kn = {'kn','clearness_dni','direct_clearness_index'};
        ALIAS.kd = {'kd','fraction_diffuse','diffuse_fraction'};

        ALIAS.Patm = {'patm','pressure','pamb'};
        ALIAS.RH = {'rh','relhum','relative_humidity'};
        ALIAS.tpw = {'tpw','total_precipitable_water'};
        ALIAS.windir = {'windir','wdir','wind_direction'};
        ALIAS.vw = {'vw','ws','windspeed','wind_speed'};
        ALIAS.Ta = {'ta','tamb','temperature_air'};

        ALIAS.albedo = {'albedo','rho'};
        ALIAS.USW = {'usw','upsw','upshw','upwelling'};
        ALIAS.soiling = {'soiling'};

        ALIAS.clearsky = {'clearsky','clear_sky','clear_samples'};
        ALIAS.CSGHI = {'csghi','clearsky_ghi','clear_sky_ghi'};
        ALIAS.CSBNI = {'csbni','clearsky_bni','clear_sky_bni','csdni','clearsky_dni','clear_sky_dni'};
        ALIAS.CSDHI= {'csdh','clearsky_dhi','clear_sky_dhi'};

        ALIAS.TL = {'TL','linke_turbidity'};
        ALIAS.AMa = {'AMa','airmass_absolute','am'};
        ALIAS.AOD = {'aod'};

        ALIAS.Ts = {'Ts','sensor_temperature'};
        ALIAS.Tm = {'Tm','module_temperature'};

        [alias,varargin] = getpairedoptions(varargin,ALIAS);

        if isa(S,'MeteoData')
            [S.data,Idx] = renamefields(S.data,alias,'-ignorecase','cat',2,varargin{:});
            [varargout{1:2}] = deal(S,Idx);
        elseif isstruct(S) || isobject(S)
            [varargout{1:2}] = renamefields(S,alias,'-ignorecase','cat',2,varargin{:});
        else
            [varargout{1:2}] = parselist(S,alias,varargin{:});
        end
    end
end
