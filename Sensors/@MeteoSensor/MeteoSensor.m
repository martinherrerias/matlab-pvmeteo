classdef MeteoSensor
% The pupose of this class is twofold:
%
%   1. MeteoSensor objects can be used by the MeteoData constructor to map data headers/fields
%      (sensor_obj.ID) into standardized variable names (sensor_obj.type).
%   2. Objects contain meta-data for the specific sensor, required to determine the sensor's
%       orientation and uncertainty (possibly as a function of solar position and weather):
%
%     obj.location - can be used to get a sensor's absolute position given relative coordinates
%     obj.mount - A structure with fields tilt, azimuth, [tracking type], etc. to be read by
%       mountrotations, to define the sensors orientation [as a function of solar position].
%
%     obj.zero_offset - [W/m²] scale-independent uncertainty, typically 3-15 W/m²
%     obj.calibration - [percent] nonlinearity, stability, spectral response, etc. Typ. 3-5 %
%     obj.temp_response - [percent/K] temperature response, typ. 0.015 - 0.25. ~ | Ta - 25°C |
%     obj.azimuth_response - [percent] directional response, proportional to BNI·cos(AOI)
%     obj.cosine_response - [W/m²] for 1000 W/m² @ 80° AOI. Typ. 5-20 W/m², ~ BNI·sin(AOI)
%     obj.tilt_response - [percent] added uncertainty for tilted sensors. Typ. 0.5 - 2.0 %
%
% The class is designed for the following workflow:
%
%   a. Use MeteoSensor.readsensorspecs to retrieve sensor information from a site-configuration
%      file ('site.sensors'). This will involve calls to the class constructor and to the lookup
%      function MeteoSensor.sensor_defaults which reads from 'sensor_defaults.xlsx'
%
%   b. Pass the resulting objects to MeteoData, to aid in data parsing and quality-control.
%
%   c. Calls to MeteoSensor.uncertainty (from within MeteoData) for advanced quality checks and
%      sensor fusion (Optimal State Estimation).
%
% The uncertainty model borrows from references [1]-[9], although many of the default values in 
% 'sensor_defaults.xlsx' had to be based on comparison with other sensors, and own judgement.
%
% References:
% [1] M. Sengupta, A. Habte, C. Gueymard, S. Wilbert, and D. Renné, “Best Practices Handbook for
%    the Collection and Use of Solar Resource Data for Solar Energy Applications: Second Edition,” 
%    NREL, NREL/TP-5D00-68886, Dec. 2017.
% [2] Driesse, A., Zaaiman, W., Riley, D., Taylor, N., Stein, J.S., 2015. Indoor and Outdoor 
%    Evaluation of Global Irradiance Sensors. Presented at the 31st European Photovoltaic Solar 
%    Energy Conference, Hamburg, Germany, p. 6.
% [3] NREL BMS Radiometer Data Uncertainty Analysis Worksheet, retrieved: 08/20/2014,
%     https://midcdmz.nrel.gov/apps/sitehome.pl?site=BMS
% [4] World Meteorological Organization. 2008. Guide to Meteorological Instruments and Methods 
%    of Observation. Geneva, Switzerland: World Meteorological Organization.
% [5] D. L. King and D. R. Myers, “Silicon-photodiode pyranometers: operational characteristics,
%   historical experiences, and new calibration procedures,” in Conference Record of the Twenty
%   Sixth IEEE Photovoltaic Specialists Conference - 1997, Anaheim, CA, USA, 1997, pp. 1285–1288.
%   doi: 10.1109/PVSC.1997.654323.
% [6] F. Vignola, “Removing Systematic Errors from Rotating Shadowband Pyranometer Data.” 2015. 
%   Accessed: Aug. 10, 2021. URL: https://www.researchgate.net/publication/288596651
% [7] S. Wilbert et al., “Uncertainty of rotating shadowband irradiometers and Si-pyranometers 
%   including the spectral irradiance error,” Cape Town, South Africa, 2016, p. 150009. 
%   doi: 10.1063/1.4949241.
% [8] W. Marion et al., “User’s Manual for Data for Validating Models for PV Module Performance,” 
%     NREL/TP-5200-61610, 1130632, Apr. 2014.
% [9] Kipp & Zonen CMP/CMA series manual, 2006
%
% See also: MeteoData, MeteoQC, MeteoSensor.readsensorspecs, MeteoSensor.uncertainty 
    
properties
    ID (1,:) char             % unique ID (to match with meteo-data header)
    type (1,:) char           % { GHI / BNI / DHI / USW }
    model (1,:) char          % custom or char-key for defaults, e.g. 'MS80','RSI',..'A','B','S','P'
    group (1,:) char          % unique ID (to match with other, related sensors)
    info (:,1) cell           % Free text field
    % units (1,:) char ?

    location = struct.empty;      % location of the instrument (absolute coords. or module coords.)
    mount struct = struct.empty;  % tilt, azimuth, type... should work with MOUNTROTATIONS

    % These are used by MeteoSensor.uncertainty
    zero_offset double = []; 
    calibration double = [];
    temp_response double = [];
    azimuth_response double = [];
    cosine_response double = [];
    tilt_response double = [];

        % % Only in new MATLAB versions...
        % zero_offset double {mustBeScalarOrEmpty,mustBeNonnegative} = []; 
        % calibration double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
        % temp_response double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
        % azimuth_response double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
        % cosine_response double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
        % tilt_response double {mustBeScalarOrEmpty,mustBeNonnegative} = [];
    
    fIAM  % function handle, used for IAM correction (see checkIAM)
    % shading ShadingResults {mustBeScalarOrEmpty} TODO
end
properties (Constant = true, Hidden = true)
    % known types and their aliases
    TYPES = {'GHI','GHI';'DHI','DHI';'BNI',{'DNI','BNI'};'GTI','GTI';...
             'USW',{'USW','UPSW','UPSHW','UP','albedo','upwelling'}};
    
    % fields recognized by constructor, that are put into mount property
    MOUNT_FLD = {'mount_type','tilt','azimuth','backtracking',...
                     'axisoffset','groundcoverratio','tracklimits'}';
                 
    AZCONV = 'eq0'; % internal (default) azimuth convention
end
methods (Hidden = true, Static = true)
    function [typ,msg] = checktype(typ)
        try
            [typ,ix] = parselist(typ,MeteoSensor.TYPES,'-soft');
            if ~all(ix > 0)
                typ = MeteoData.findfields(typ);
            end
            msg = '';
        catch ERR
            if ~strcmp(ERR.identifier,'parselist:unknown'), rethrow(ERR); end
            msg = ['Unknown sensor type: ',typ];
            if nargout < 2, warning('MeteoSensor:unknown',msg); end
        end
    end
    
    function [S,msg] = sensor_defaults(type,mdl)
    % S = SENSOR_DEFAULTS(TYPE,MDL) - search in worksheet "sensor_defaults.xls" for default uncertainty
    %   parameters of a given sensor TYPE and MDL. 

        persistent DEF;
        if isempty(DEF)
            [data,txt] = xlsread('sensor_defaults.xlsx',1);
            [n,m] = size(data);
            DEF = [array2table(data,'VariableNames',txt(1,2:m+1)),...
                   cell2table(txt(2:n+1,[1,m+2:end]),'VariableNames',txt(1,[1,m+2:end]))];
            DEF.alias = cellfun(@(s) ['^(' s ')$'],DEF.alias,'unif',0);
        end

        switch type
        case 'BNI'
            D = DEF(DEF.B ~= 0,:);
        case {'GHI','DHI'}
            a = type(1);
            D = DEF(DEF.(a) == 1,:);
            D.tilt_response(:) = 0;
        case {'GTI'}
            D = DEF(DEF.G == 1,:);
        case {'USW'}
            D = DEF(DEF.D == 1,:);
        otherwise
            S = struct();
            msg = ['No default properties for sensor type: ',type];
            if nargout < 2, warning(msg); end
            return;
        end
        assert(~isempty(D),'Failed to read sensor defaults');

        try
            assert(~isempty(mdl),'Unknown sensor model');
            [~,idx] = parselist(mdl,[D.model,D.alias],'sensor model','-regexp','-soft');
            if numel(idx) > 1
                idx = find(idx,1);
                msg = sprintf('Found more than one match for "%s", using %s',mdl,D.model{idx});
            elseif idx == 0
                error('Model "%s" not found',mdl);
            else
                msg = '';
            end
        catch ERR
            idx = contains(D.type_default,type);
            if nnz(idx) == 1
                msg = sprintf('%s: using default (%s) for %s',ERR.message,D.model{idx},type);
            else
                S = struct();
                msg = ['No default model for sensor type: ',type];
                if nargout < 2, warning(msg); end
                return;
            end
        end

        S = rmfield(table2struct(D(idx,:)),{'alias','G','D','B','type_default'});
        if ~isempty(S.fIAM), S.fIAM = str2func(S.fIAM); end
        
        if nargout < 2, warning(msg); end
    end
end
methods
    function [obj,msg] = MeteoSensor(varargin)
    % S = MeteoSensor(X) - parse structure/object X into a MeteoSensor. Complete missing fields
    %   with defaults, whenever available.
    % S = MeteoSensor(KEY) - searches character string KEY for recognized sensor types, and takes
    %   the rest as model. e.g. 'GHI_PSP' is equivalent to struct('type','GHI','model','PSP')
    % S = MeteoSensor(..,'name',val) - provide/modify individual properties
    %
    % See also: MeteoSensor.readsensorspecs
           
        if nargin == 0, return; end
        fields = [fieldnames(obj); MeteoSensor.MOUNT_FLD];
                
        [S,varargin] = getpairedoptions(varargin,fields);
        if ~isempty(varargin)
            assert(numel(varargin) == 1,'Unrecognized argument(s)');
            
            if isstruct(varargin{1}) || isobject(varargin{1})
            % IRRADIANCESENSOR(S,..)
                if istable(varargin{1}), varargin{1} = table2struct(varargin{1}); end
                S = completestruct(S,varargin{1});
                
            elseif ischar(varargin{1})
            % IRRADIANCESENSOR(KEY,..), for KEY = 'TYPE_MODEL'
                key = varargin{1};
                alltypes = cat(2,obj.TYPES{:,2});
                [typ,mdl] = regexpi(key,strjoin(alltypes,'|'),'match','split');
                mdl = regexprep([mdl{:}],'(^[\W_]*)|([\W_]*$)','');
                if ~isempty(typ), S.type = typ{:}; end
                if ~isempty(mdl), S.model = mdl; end
            else
                error('Expecting structure, MeteoSensor, or character key');
            end
        end
        if isfield(S,'info'), S.info = cellstr(S.info); end
        
        % Parse mounting information
        mountfld = MeteoSensor.MOUNT_FLD;
        mountfld(cellfun(@(f) ~isfield(S,f) || isempty(S.(f)),mountfld)) = [];
        if ~isempty(mountfld)
            for f = mountfld'
                S.mount(1).(strrep(f{1},'mount_','')) = S.(f{1}); 
            end
            S = rmfield(S,mountfld);
            if isfield(S.mount,'azimuth')
                if ischar(S.mount.azimuth)
                    c = regexpi(S.mount.azimuth,'(.*)(([NSEW]2[NSEW])|(eq[+-\d.]*)).*','tokens');
                    if ~isempty(c)
                        S.mount.azimuth = str2double(c{1}{1});
                        conv = c{1}{2};
                        S.mount.azimuth = solarposition.fixazimuth(S.mount.azimuth,conv,'N2E');
                    else
                        S.mount.azimuth = str2double(S.mount.azimuth);
                    end
                end
            end
        end
        if ismember(S.type,{'GHI','DHI','BNI','GTI','USW'})
            if ~isfield(S,'mount') ||  isempty(S.mount), S.mount = struct(); end
        end
        switch S.type
        case 'BNI'
            S.mount = completestruct(S.mount,struct('type','2a'));
        case {'DHI','GHI'}
            S.mount = completestruct(S.mount,struct('type','0a','tilt',0,'azimuth',0,'slope',0));
        case 'GTI'
            S.mount = completestruct(S.mount,struct('type','0a','tilt',NaN,'azimuth',NaN,'slope',0));
        case 'USW'
            S.mount = completestruct(S.mount,struct('type','0a','tilt',180,'azimuth',0,'slope',0));
        end

        PARAMS = {'zero_offset','calibration','temp_response','azimuth_response',...
                   'cosine_response','tilt_response'};
        
        missing = true(size(fields));
        for j = 1:numel(fields)
            if ~isfield(S,fields{j}) || isempty(S.(fields{j})), continue; end
            obj.(fields{j}) = S.(fields{j});
            if contains(fields{j},PARAMS)
                validateattributes(obj.(fields{j}),{'numeric'},{'scalar','nonnegative'},'',fields{j})
            end
            missing(j) = false;
        end

        [obj.type,msg] = MeteoSensor.checktype(obj.type);
        if ~isempty(msg)
            obj.info{end+1} = msg;
            if nargout < 2, warning('For sensor %s: %s',obj.ID,msg); end
            return;
        end

        if ~any(missing), return; end
        
        % Complete with default parameter values based on type and model
        if ~isfield(S,'model'), S.model = ''; end
        lastwarn('');
        [DEF,msg] = MeteoSensor.sensor_defaults(S.type,S.model);
        if ~isempty(msg)
            obj.info{end+1} = msg;
            if nargout < 2
                if ~isempty(obj.ID)
                   msg = sprintf('%s sensor %s',msg,obj.ID);
                end
                warning(msg); 
            end
        end

        missing = missing & isfield(DEF,fields);
        if any(missing)
            obj.info{end+1} = sprintf('Using %s from default "%s"',...
                shortliststr(fields(missing),'field',Inf,'quotes',''''),DEF.model);
            for j = find(missing)'
                obj.(fields{j}) = DEF.(fields{j}); 
            end
        end
    end

    U = uncertainty(S,x,AOI,B,Ts)
    
    function S = parselocations(S,AbsLoc,Trck)
    % S = MeteoSensor.parselocations(S,Loc,Trck) 
    % Parse relative sensor location labels of the type S(j).location = 'N + [dx,dy,dz]', where
    % N is a mount ID matching TRCK and [dx,dy,dz] an offset in mount coords. Return absolute  
    % sensor location S(j).location = [x,y,z], and angles S(j).tilt, S(j).azimuth.
        
        if nargin < 3, Trck = []; end
    
        loc = {S.location};
        relative = cellfun(@ischar,loc);
        known = cellfun(@(x) isnumeric(x) && any(numel(x) == [2,3]),loc);

        if ~isempty(AbsLoc) && ~(isstruct(AbsLoc) && isempty(fieldnames(AbsLoc)))
            try
                [AbsLoc.longitude,AbsLoc.latitude,AbsLoc.altitude] = ...
                    checkcoordsystem(0,0,0,AbsLoc,'output','abs','-inland');
                AbsLoc = parselocation(AbsLoc);
            catch
                warning('Invalid location/coordinate system')
                AbsLoc = [];
            end
        else
           AbsLoc = []; 
        end
        
        if any(known)
            if ~isempty(Trck)
                [~,~,~,TLoc] = checkcoordsystem(0,0,0,Trck,'output','prj');
            end
            d = solarposition.arcdist(AbsLoc.latitude,AbsLoc.longitude,...
                                      TLoc.origin(2),TLoc.origin(1),6370);
            if d > 1
                warning(['Sensor origin seems %0.1f km away from Project, ',...
                         'resulting project coordinates might be meaningless'],d);
            end
        end
        
        if ~isempty(AbsLoc) && any(known)
            for j = find(known)
                [x,y,z] = checkcoordsystem(loc{j}(1),loc{j}(2),loc{j}(end:3),AbsLoc,'output','abs');
                if ~isempty(Trck)
                    [x,y] = abs2prj(x,y,TLoc.origin(2),TLoc.origin(1),TLoc.rotation);
                    z = z - AbsLoc.altitude + TLoc.altitude;
                end
                S(j).location = [x,y,z];
            end
        end

        % parse relative locations: 'n + [dx,dy,dz]'
        L = regexp(loc(relative),'(?<trck>\d+)[\s+]{1,3}(?<offset>\[.+\])','names');
        relative(relative) = ~cellfun(@isempty,L);

        if any(relative)
            try
                assert(~isempty(Trck),'Missing mounts structure');
                mountrotations(Trck,0,0,'N2E');
            catch ERR
                warning('Invalid layout structure, failed to parse relative locations: \n%s',...
                    getReport(ERR,'basic'));
                relative(:) = false;
            end
        end
        if any(relative)
            FLD = [MeteoSensor.MOUNT_FLD;{'type';'centers';'slope'}];
            Trck = rmfield(Trck,setdiff(fieldnames(Trck),FLD));
            if ~isfield(Trck,'axisoffset'), Trck.axisoffset = 0; end
            
            Nu = size(Trck.centers,2);
            sidx = find(relative);
            L = cat(1,L{:});
            for j = 1:numel(L)
                k = sidx(j);
                offset = str2num(L(j).offset); %#ok<ST2NM>
                offset(end:3) = 0;
                tidx = str2double(L(j).trck);
            
                Mt = filterstructure(Trck,tidx,Nu);
                Mt = filterstructure(Mt,tidx,Nu,'dim',2);
                Mt.axisoffset = offset(:) + Mt.axisoffset(:);
            
                if strcmp(Mt.type,'0a')
                    [R,ZXZ] = mountrotations(Mt,0,0,'N2E');
                    S(k).location = Mt.centers + R*Mt.axisoffset;
                    S(k).mount(1).tilt = ZXZ(1,1,2);
                    S(k).mount(1).azimuth = -ZXZ(1,1,1); % 'N2W' -> 'N2E'
                    S(k).mount(1).slope = 0;
                    S(k).mount(1).type = '0a';
                else
                    S(k).location = Mt.centers;
                    S(k).mount = rmfield(Mt,'centers');
                end
            end
        end
    end
    
    function S = checkIAM(S)
        warning_resetter = naptime('completestruct:BnotinA','error'); %#ok<NASGU>
        
        msg = cell(1,numel(S));
        for j = 1:numel(S)
            if ~isempty(S(j).fIAM)
                try 
                    [S(j).fIAM,~,~,msg{j}] = checkIAM(S(j).fIAM);
                    continue;
                catch ERR
                    msg{j} = ['Invalid IAM for sensor ''' S(j).ID ''': ' ERR.message];
                    S(j).fIAM = []; 
                end
            else
                msg{j} = ['No IAM for sensor ''' S(j).ID ''''];
            end

            switch S(j).type
            case 'BNI'
                S(j).fIAM = checkIAM('hcpv','MaxIncAngle',2.5);
                msg{j} = [msg{j} ', assuming ideal 2.5° aperture angle'];
          
            case {'GTI','GHI','USW','DHI'} % assume ideal response for DHI/USW?
                
                if S(j).cosine_response > 0
                    try
                        [S(j).fIAM,txt] = checkIAM('martinruiz','maxloss',S(j).cosine_response/1000);
                        msg{j} = [msg{j} ', estimating' txt 'from cosine-response'];
                    end
                end
                if isempty(S(j).fIAM)
                    S(j).fIAM = checkIAM('none');
                    msg{j} = [msg{j} ', assuming ideal response'];
                end
                
            otherwise, msg{j} = '';
            end 
        end
        msg = strjoin(msg(~cellfun(@isempty,msg)),newline());
        if ~isempty(msg)
            warning(msg)
        end
    end
    
    function R = struct(S)
        R = cell2struct(cell(S),fieldnames(S));
    end
    function C = cell(S)
        fld = fieldnames(S);
        C = arrayfun(@(s) cellfun(@(f) {s.(f)},fld),S,'unif',0);
        C = cat(2,C{:});
    end
    function T = table(S)
        T = cell2table(cell(S)','VariableNames',fieldnames(S));
    end
end
methods (Static = true)
    [S,Loc,U] = readsensorspecs(specsfile)
end
end
