function [S,Loc,U] = readsensorspecs(specsfile)
% [S,LOC,U] = READSENSORSPECS(YMLFILE) - Read a YAML file with sensor specifications, and
%   parse output into an array of MeteoSensor objects.
%
%   Basic block definition (see [1][2] and 'example.sensors' for format details):
%
%     - # Comment/sensor-description (space-dash + indendation mark new sensor)
%       # CAUTION: tabs (\t) for indentation are forbidden in YAML! replace with spaces
%
%         id: GHI_1         # unique valid variable name, matching *.meteo column header
%         type: GHI         # See recognized keys in MeteoSensor.uncertainty
%         model: S          # See recognized defaults in "sensor_defaults.xls"
%         zero_offset: 5.0  # ...
%         ... 
%         tilt_response: 0.0
%         location: [0,0,2]  # static project coordinates [x,y,z], or relative position ...
%                            # N + [dx,dy,dz], where N is a mount ID (defined in *.mounts file),
%                            # and [dx,dy,dz] an offset in mount coords.
%         tilt: 20           # Degrees from horizontal (not required for relative position).
%         azimuth: 10.0 N2E  # Degrees [convention]
%
%         members:          # sensor/group members inherit all properties from the parent,
%             -
%             id: GHI_1_sensor_temperature
%             type: Ts
% 
% NOTES:
%   
%   Name YMLFILE with *.sensors extension and place on PWD for automatic detection.
%
%   Each sensor must contain a sensor TYPE and a unique ID. Structures that don't match these
%   criteria will cause a warning, and can be retrieved as third output argument U.
%
%   Sensor types that are not recognized (see MeteoData.varnames) will cause a warning, but
%   be kept. They can be used to cluster variables in MeteoData, as long as there are other 
%   sensors with the same (case sensitive) type.
%
%   Location information will be used to get sensor shading factors (TODO). Right now they
%   can be used to get the orientation of field-mounted sensors (MeteoSensor.parselocations). 
%
%   Tilt and azimuth are important for static GTI sensors, can be left out otherwise.
%
% REF:
%   [1] "YAMLMatlab". Jiri Cigler et al. 2012. https://code.google.com/archive/p/yamlmatlab/
%   [2] "YAML Ain’t Markup Language (YAML™) Version 1.2", https://yaml.org/spec/1.2/spec.html
%
% See also: MeteoSensor, MeteoData, MeteoSensor.uncertainty

    import yaml.*
    PROPS = fieldnames(MeteoSensor());
    PROPS = [PROPS; 'members'; 'tilt'; 'azimuth'];
    % PROPS = [PROPS; cellfun(@(p) ['members.' p],PROPS,'unif',0)];

    notemptystr = @(x) ischar(x) || isstring(x) && ~isempty(x);
    grouplist = {};

    report = {};

    specsfile = pickfile(specsfile);
    S = ReadYaml(specsfile,1,true);

    if isscalar(S) && isstruct(S)
        S = renamefields(S,{'location','sensors'}');
        if isfield(S,'location')
            try Loc = parselocation(S.location);
            catch, Loc = [];
            end
        end
        if isfield(S,'sensors'), S = S.sensors; end
    else
        Loc = [];
    end
    [S,valid] = parsesensors(S);

    U = S(~valid);
    S = cat(1,S{valid});

    if ~all(valid) && nargout < 3
       report{end+1} = sprintf('Ignoring %d invalid/incomplete maps: \n%s',nnz(~valid),...
            strjoin(cellfun(@dispnested,U,'unif',0),'\n\n'));
    end              
    if ~isempty(report)
        report = uniquecell(report,'stable');
        warning(strjoin(report,'\n')); 
    end

    function [S,valid] = parsesensors(X,parentgroup)
    % Parse YAML output into MeteoSensor objects, copy parent properties into members, 
    % and replace hierarchy with group tags.

        if nargin < 2, parentgroup = ''; end
        X = X(:);

        if iscell(X)
            [S,valid] = cellfun(@(x) parsesensors(x,parentgroup),X,'unif',0);
            S = cat(1,S{:});
            valid = cat(1,valid{:});
            return;
        end

        S = cell(numel(X),1);
        valid = cell(numel(X),1);

        for j = 1:numel(X)
            s = X(j);
            assert(isstruct(s),'Expecting structure or cell-array of structures');

            [s,~,R] = renamefields(s,PROPS,'-ignorecase');
            good = all(isfield(s,{'ID','type'})) && notemptystr(s.ID) && notemptystr(s.type);

            % Make sure there is a 'group' field 
            if ~(isfield(s,'group') && notemptystr(s.group))
                if ~isempty(parentgroup)
                    s.group = parentgroup;
                elseif isfield(s,'ID') && notemptystr(s.ID)
                    s.group = s.ID;
                else
                    s.group = matlab.lang.makeUniqueStrings('unnamed_group',grouplist);
                    grouplist{end+1} = s.group; %#ok<AGROW>
                end
            elseif ~isempty(parentgroup)
            % ... and inherit parent's group(s), if available

                if ~contains(parentgroup,s.group,'ignorecase',true)
                    s.group = [s.group '.' parentgroup];
                end
            end

            if isfield(s,'members')
            % Recursive parse of children, inheriting anything missing

                fill = @(x) completestruct(renamefields(x,PROPS,'-ignorecase'),...
                    rmfield(s,{'members','group'}));
                if isstruct(s.members)
                    s.members = arrayfun(fill,s.members,'unif',0);
                elseif iscell(s.members)
                    s.members = cellfun(fill,s.members,'unif',0);
                end
                if isfield(s,'ID') && notemptystr(s.ID) && ...
                        ~contains(s.group,s.ID,'ignorecase',true)
                    grp = [s.group '.' s.ID];
                else
                    grp = s.group; 
                end

                [s.members,goodchildren] = parsesensors(s.members,grp);
            else
                goodchildren = false(0,0); 
            end

            if ~good && ~any(goodchildren), R = completestruct(R,s); end

            if ~isempty(R) && ~isempty(fieldnames(R))
                msg = sprintf('unrecognized fields: \n%s',dispnested(R));
                if good
                    if ~isfield(s,'info'), s.info = {}; end
                    s.info = cellstr(s.info);
                    s.info{end+1} = msg;
                    warning('In %s: %s',s.ID,msg);
                else
                    if isfield(s,'group') && notemptystr(s.group)
                        msg = sprintf('In %s: %s',s.group,msg);
                    else
                        msg = strrep(msg,'unrecognized fields:','Incomplete structure:');
                    end
                    warning(msg);
                end
            end

            % destroy family hierarchy - down with the patriarchy!
            if good
                if any(goodchildren)
                    if ~iscell(s.members),s.members = num2cell(s.members); end
                    S{j} = [{s};s.members];
                    valid{j} = [true;goodchildren];
                else
                    S{j} = {s};
                    valid{j} = true;
                end
            else
                if any(goodchildren)
                    if ~iscell(s.members),s.members = num2cell(s.members); end
                    S{j} = s.members(goodchildren);
                    valid{j} = goodchildren;
                else
                    S{j} = {R};
                    valid{j} = false;
                end                    
            end                 
        end

        valid = cat(1,valid{:});
        S = cat(1,S{:}); 

        for j = 1:numel(valid)
            if ~valid(j) || ~isstruct(S{j}), continue; end
            try
                [S{j},msg] = MeteoSensor(S{j});
                if ~isempty(msg)
                    report{end+1} = msg; %#ok<AGROW>
                end
            catch ERR
                valid(j) = false;
                report{end+1} = sprintf('In %s: %s',S{j}.ID,ERR.message); %#ok<AGROW>
            end
        end  

    end
end 