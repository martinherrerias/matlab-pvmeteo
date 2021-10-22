 function MD = checksources(MD,varargin)
% MD = checksources(MD) - verify that MD.source is a cell array of cellstrings, one for
%   each column, i.e. size(MD(:,j),2) = numel(MD.source{j}). Fill with defaults if missing.
%   For a variable A with size(MD.A,2) == m, default is {'A.1','A.2',..,'A.m'}.
%
%   Check that sources are unique (they serve as ID's for sensors/models).
%
% MD = checksources(MD,IDX1,IDX2,..) - combine rename indices IDX1, IDX2 (structures with
%   form IDX.target = source as returned by RENAMEFIELDS), into MD.source:
%
%       Before: fieldnames(MD) -- ... > IDX2 --> IDX1 --> original IDs
%       After:  fieldnames(MD) --> MD.source (original IDs)

    MD.data.Properties.DimensionNames{1} = 't';

    varnames = MD.data.Properties.VariableNames';
    n = cellfun(@(f) size(MD.data.(f),2),varnames);

    default = @(c,n) arrayfun(@(k) sprintf('%s.%d',c,k),1:n,'unif',0);

    if ~isprop(MD.data.Properties.CustomProperties,'source')
       MD.data = addprop(MD.data,'source','variable');
    end

    if ~isempty(MD.source)
    % Use existing aliases as an additional layer, to be merged with defaults.
        validateattributes(MD.source,{'cell'},{'vector','numel',numel(varnames)});
        varargin{end+1} = cell2struct(MD.source',varnames);
    end

    % ... then reset to default (trivial) sources: A: {A.1,A.2,..}, ..
    MD.data.Properties.CustomProperties.source = cellfun(default,varnames,num2cell(n),'unif',0);

    source_list = MD.source;
    source_list(n == 1) = cellfun(@cellstr,source_list(n == 1),'unif',0);
    source_list = cat(2,source_list{:}); 

    [~,ic] = unique(source_list);
    if numel(ic) ~= numel(source_list)
        ic = setdiff(1:numel(source_list),ic); % repeated
        error('Non-unique source %s',shortliststr(source_list(ic),'ID','colon',':'));
    end

    if numel(varargin) == 0, return; end

    for k = numel(varargin):-1:1

        Idx = varargin{k};            
        try
            validateattributes(Idx,{'struct'},{'scalar'},'',sprintf('Idx{%d}',k));
            fld = fieldnames(Idx);
            Idx = struct2cell(Idx);
        catch
            error('Invalid alias, expecting structure with char/cellstr fields');
        end

        timefields = strcmp(fld,'t');
        fld(timefields) = [];
        Idx(timefields) = [];

        m = cellfun(@(x) 1 + iscell(x)*(numel(x)-1),Idx,'unif',0);
        trivial = cellfun(@(f,a,n) isempty(a) || isequal(f,a) || ...
                                                 isequal(default(f,n),a),fld,Idx,m);
        if all(trivial), continue; end
        fld(trivial) = [];
        Idx(trivial) = [];

        for j = 1:numel(fld)
            if iscell(Idx{j})
                fld{j} = default(fld{j},numel(Idx{j}));
            else
                fld{j} = cellstr(fld{j});
                Idx{j} = cellstr(Idx{j});
            end
            if isscalar(fld{j}) && isfield(MD,fld{j}) && ...
               numel(cellstr(Idx{j})) == size(MD.data.(fld{j}{1}),2)    % simple rename
                ib = sum(n(1:fieldidx(MD,fld{j})));
                ix = true(1,numel(Idx{j}));
            else
                [ix,ib] = ismember(fld{j},source_list);
            end
            if ~all(ix)
                warning('Index %s not found in existing sources',shortliststr(fld{j},'field'));
            end
            source_list(ib(ix)) = Idx{j}(ix);
        end
    end

    MD.data.Properties.CustomProperties.source = mat2cell(source_list,1,n);
end