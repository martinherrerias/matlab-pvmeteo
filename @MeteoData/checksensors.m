function MD = checksensors(MD,sensors,assumedefaults)
% MD = CHECKSENSORS(MD,[S,ASSUMEDEFAULTS]) - Parse existing MD.sensors, and optionally add new
%   sensors S. Both MD.sensors and S are expected to be vectors of METEOSENSOR objects, with
%   unique IDs
%
% See: METEOSENSOR, METEODATA, METEODATA.CHECKSOURCES, METEODATA.UNCERTAINTY

    % MD.data.Properties.CustomProperties.sensors = sensors; 

    if nargin < 3, assumedefaults = false; end

    varnames = MD.data.Properties.VariableNames;
    n = cellfun(@(f) size(MD.data.(f),2),varnames);

    if nargin < 2 || isempty(sensors), sensors = MeteoSensor.empty; end
    if isempty(MD.sensors), MD.sensors = MeteoSensor.empty; end
    
    validateattributes(sensors,{'MeteoSensor'},{},'','sensors');
    validateattributes(MD.sensors,{'MeteoSensor'},{},'','MD.sensors');

    try
        assert(all(cellfun(@numel,MD.source) == n));
        source_list = cat(2,MD.source{:}); 
        assert(iscellstr(source_list) && numel(unique(source_list)) == numel(source_list)); %#ok<ISCLSTR>
    catch
        MD = checksources(MD);
        source_list = cat(2,MD.source{:}); 
    end

    sensIDs = {MD.sensors.ID};
    newIDs = {sensors.ID};

    % Replace existing sensors with the same ID, add new ones to list
    [iv,ib] = ismember(newIDs,sensIDs);
    if any(iv)
        warning('Replacing %s with %s',shortliststr(sensIDs(ib(iv)),'sensor'),...
                                       shortliststr(newIDs(iv),'sensor'));
        MD.sensors(ib(iv)) = sensors(iv);
    end
    MD.sensors = cat(2,MD.sensors,sensors(~iv));
    sensIDs = {MD.sensors.ID};

    [~,ic] = unique(sensIDs);
    if numel(ic) ~= numel(sensIDs)
        ic = setdiff(1:numel(sensIDs),ic); % repeated
        error('Non-unique sensor %s',shortliststr(sensIDs(ic),'ID','colon',':'));
    end

    [assigned_sensors,ib] = ismember(sensIDs,source_list);
    assigned_sources = false(size(source_list));
    assigned_sources(ib(assigned_sensors)) = true;

    % Use unmatched sensors ID and type as dictionary to rename variables
    if ~all(assigned_sensors)
        assigned_vars = cellfun(@any,mat2cell(assigned_sources,1,n));
        % known = ismember(varnames,MeteoData.varnames('all'));
        
        [~,ix] = parselist(sensIDs(~assigned_sensors)',varnames(~assigned_vars)','-soft');
        if iscell(ix)
            [~,ix] = parselist(sensIDs(~assigned_sensors),varnames(~assigned_vars),'-soft','-matchcase');
        end
        if any(ix > 0)
            
            
            is = find(~assigned_sensors); is = is(ix > 0);
            iv = find(~assigned_vars); iv = iv(ix(ix > 0));
            
            newTypes = {MD.sensors(is).type};
            newIDs = varnames(iv);
            
            silly = (n(iv) == 1) & strcmpi(newTypes,newIDs); % source = scalar variable name
            if any(silly)
                newIDs(silly) = cellfun(@(x) [x,'.1'],newIDs(silly),'unif',0);
                newIDs(silly) = matlab.lang.makeUniqueStrings(newIDs(silly),cat(2,MD.source{assigned_vars}));
                warning('Changing to sensor %s to avoid confusion with variable %s',...
                    shortliststr(newIDs(silly),'ID'),shortliststr(newTypes(silly),'name'));
                [MD.sensors(is(silly)).ID] = deal(newIDs{silly});
                MD.data.Properties.CustomProperties.source(iv) = newIDs(silly);
                sensIDs = {MD.sensors.ID};
            end
            if any(~silly)
                [MD.data,Idx] = renamefields(MD.data,[newIDs(~silly);newTypes(~silly)]','cat',2);
                MD = checksources(MD,Idx);
            else
                MD = checksources(MD);
            end

            [MD.sensors(is).ID] = deal(newIDs{:});

            source_list = cat(2,MD.source{:}); 
            varnames = MD.data.Properties.VariableNames;
            n = cellfun(@(f) size(MD.data.(f),2),varnames);

            [assigned_sensors,ib] = ismember(sensIDs,source_list);
            assigned_sources = false(size(source_list));
            assigned_sources(ib(assigned_sensors)) = true;
        end
    end

    % % attempt to deal unassigned sensors by type
    % if ~all(ia)
    %     newTypes = {MD.sensors(~ia).type};
    %     repelem(varnames,n)
    %     ...
    % end

    if any(~assigned_sensors)
        warning(shortliststr(sensIDs(~assigned_sensors),'Removing unassigned sensor'));
        MD.sensors(~assigned_sensors) = [];
    end
    
    [areflags,parents] = findflagsensors(MD.sensors,source_list);
    if any(areflags)
        [MD,X] = rmsource(MD,{MD.sensors(areflags).ID});
        MD = parseflags(MD,X,parents);
        MD.sensors(areflags) = [];
    end

    knowntype = ismember(varnames,MeteoSensor.TYPES(:,1));
    assumed = repelem(knowntype,n) & ~assigned_sources;

    if ~any(assumed), return; end

    if ~assumedefaults
        warning('checksensors:missing',shortliststr(source_list(assumed),'Missing sensor'));
        return;
    end

    assumed_vars = cellfun(@any,mat2cell(assumed,1,n));
    typeidx = assumed.*repelem(1:numel(varnames),n);
    
    for j = find(assumed_vars)
        idx = (typeidx == j);
        S = cellfun(@(id) MeteoSensor(varnames{j},'ID',id),source_list(idx));
        [S.ID] = deal(source_list{idx});
        MD.sensors = [MD.sensors(:)',S];
    end
end

function [areflags,parents] = findflagsensors(sensors,source_list)

    assigned = ismember({sensors.ID},source_list);

    areflags = contains({sensors.type},'flag');
    areflags = areflags & assigned;
    if ~any(areflags), parents = {}; return; end

    parents = {sensors(areflags).group};
    [found,ib] = ismember(parents,source_list);

    if all(found), parents = source_list(ib); return; end
    
    for j = 1:numel(source_list)
        match = endsWith(parents(~found),['.' source_list{j}]);
        if nnz(match) == 1
            ib(~found) = match.*j;
            found(~found) = match;
            if all(found), break; end
        end
    end
    
    areflags(~found) = false;    
    parents = source_list(ib(found));
end

function MD = parseflags(MD,X,parents)

    MAXBITS = meteoQC.BITS - numel(meteoQC.FLAGS);
    validateattributes(X,{'numeric','logical'},{'2d','integer','nonnegative','finite'},'','flags');
    
    [u,~,ia] = unique(X(:),'sorted');
    if u(1) == 0 && mode(X(:)) == 0
        u = u(2:end);
        ia = ia -1;
    end
    nu = nnz(u);
    
    if nu <= MAXBITS
    % Flag by unique IDs: label flags as 'u1','u94',... for any existing numbers 23,94,
        B = sparse(1:numel(X),ia,true);
        u = arrayfun(@(j) sprintf('u%d',j),u,'unif',0);
    else
        nu = floor(log2(max(abs(u))))+1;
        B = bitget(repmat(u,1,nbit),repmat(1:nu,numel(u),1),meteoQC.TYPE);
        used = any(B,1);
        nu = nnz(used);
        if nu > MAXBITS
            warning('No space for new flags: setting all existing flags as "other"');
            B = (X(:) > 0);
            u = {'other'};
        else
        % 
            B = B(ia,used);
            u = arrayfun(@(j) sprintf('b%d',j),find(used));
        end
    end
    
    [b,MD.flags] = flagbit(MD.flags,u);
    [~,F,iv,ic] = getbysource(MD,parents);
    
    for j = 1:numel(b)
        F(B(:,j)) = bitset(F(B(:,j)),b(j));
    end
    
    for j = unique(iv)'
        inj = (iv == j);
        fld = MD.data.Properties.VariableNames{j};
        MD.flags.(fld)(:,ic(inj)) = F(:,ic(inj));
    end
end
