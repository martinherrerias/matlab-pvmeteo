classdef MeteoData < table_ish
% Class designed to contain regularly sampled/averaged meteorological observations and related
% parameters (solar-position, clear-sky models, etc.) in the form of a TIMETABLE), along with 
% meta-data, Quality-Control flags and measurement uncertainty.
%
% Since the Tabular class is sealed (thank you, MathWorks), the class inherits from the TABLE_ISH
% wrapper. Data is stored in timetable obj.data. Variables can be accessed as MD.data.VAR, and
% outside class methods also directly as MD.VAR. Upon object contruction names are parsed to the
% internal variable naming convention (see METEODATA.WHATIS, VARNAMES, and FINDFIELDS)
%
%   TODO: switch to standard naming convention https://duramat.github.io/pv-terms/
%
% Each variable MD.VAR can be an Nt x M array of Nt measurements from M sensors of type VAR.
% Individual channels (i.e. sensors) are labeled by a unique character-ID: obj.source{j}{k} for
% the kth channel of the jth variable. See GETBYSOURCE, GETSOURCEOF, SETSOURCEOF)
%
% Channels can (and should, if possible) be tied to a METEOSENSOR object S in the obj.sensors list,
% by setting S.type == VAR, and S.ID ==  obj.source{j}{k}. This can also be used to rename columns
% automatically upon construction.
%
% Meta-data is included in class properties:
%
%     obj.t (obj.data.Properties.RowTimes) [datetime] interval/sample timestamps (see obj.interval)
%     obj.timestep (obj.data.Properties.TimeStep) scalar DURATION
%     obj.source unique column char IDs, pointing to sensor or model IDs for each variable
%
%     obj.location - structure, with fields latitude, longitude, altitude,... see PARSELOCATION
%     obj.interval - {'i','c','b','e'} interval summarization (instant / center / beginning / end)
%     
%     obj.options - structure of options used by COMPLETEMETEODATA, LOADOBJ, etc.
% 
%     obj.info - free (cell-string) field. Other information can also be swept under the rug as 
%       part of obj.data.Properties.UserData or obj.data.Properties.CustomProperties.
%
%     obj.sensors - METEOSENSOR object array
%     obj.models - FUTURE, function-handles? structure? for model uncertainty quantification
%
%     obj.flags - METEOQC object, that holds a table obj.flags.data with the same fields as 
%           obj.data with a binary-word of flags for each data point.
%
%     obj.uncertainty - table_ish object with the same structure as obj.data, with estimated 
%           expanded uncertainty for each data point with an assigned sensor.
%
% See also: METEOSENSOR, METEOQC, METEODATA.METEODATA, METEODATA.COMPLETEMETEODATA

    
properties (GetAccess = public, SetAccess = public)    
    
    location = struct.empty;  % structure, see PARSELOCATION
    interval = '';            % {'i','c','b','e'} interval summarization (instant / center / beginning / end)
    timestep                  % DURATION scalar (MD.data.Properties.TimeStep, if regular)
    info ={''}                % free cellstring
    
    options = struct.empty    % SimOptions.meteo and others
    flags = meteoQC.empty;    % MeteoQC object
    uncertainty

    sensors = MeteoSensor.empty;  % array of MeteoSensor objects
    models = {}                   % function-handles? structure?
end
properties (GetAccess = public, SetAccess = private)
    status       % structure of flags {regular,tested,complete}
end
properties (Dependent = true) % encoded in TimeTable object
    t           % MD.data.Properties.RowTimes: [datetime] interval/sample timestamps
    source      % MD.Properties.CustomProperties.source - unique column char IDs, pointing to
                % sensor or model IDs
end
properties (Dependent = true, SetAccess = private) 
    Nt
end
properties (Transient = true)
    missing
    dark
    available
end
properties (Constant, Hidden = true)
   % Fields to get (if opt.useweb) from McClear/MERRA2 files
    % VAR_NAMES = 
    VARS_MCCLEAR = {'CSGHI','CSBHI','CSDHI','CSBNI','albedo','tpw'}; % Also read, but so far not used: 'TOA','sza','fiso','fvol','fgeo'
    VARS_MERRA2 = {'Ta','RH','Patm','vw','windir'}; % Also read, but so far not used: 'rainfall','snowfall','snowdepth','GHI' 
end
methods
    function t = get.t(MD), t = MD.data.Properties.RowTimes; end
    function MD = set.t(MD,t)
        if istable(MD.data), MD.data = table2timetable(MD.data,'rowtimes',t);
        else, MD.data.Properties.RowTimes = t;
        end
    end
    
    function dt = get.timestep(MD)
        if isregular(MD.data) && MD.interval ~= 'i'
            dt = MD.data.Properties.TimeStep;
            % if ~isempty(MD.timestep) && isfinite(MD.timestep) && ~isequal(MD.timestep,dt)
            %     warning('Inconsistent timestep: %s (property) vs %s (table)',...
            %         isoduration(dt),isoduration(MD.timestep));
            % end
            MD.timestep = dt;
        end
        dt = MD.timestep;
    end
    function MD = set.timestep(MD,dt)
        if isempty(dt), dt = minutes(NaN); end
        validateattributes(dt,'duration',{'scalar'});
        MD.timestep = dt;
        if ~isempty(MD.data) && isregular(MD.data)
            dt0 = MD.data.Properties.TimeStep;
            if abs((dt - dt0)/max(dt,dt0)) > 1e-2
                warning('Inconsistent timestep: %s (property) vs %s (table)',...
                    isoduration(dt),isoduration(dt0));
            end
        end
    end
    
    function Nt = get.Nt(MD), Nt = size(MD.data,1); end
    
    function f = get.dark(MD)
        if isempty(MD.dark), MD = refresh(MD); end; f = MD.dark;
    end
    function f = get.missing(MD)
        if isempty(MD.missing), MD = refresh(MD); end; f = MD.missing;
    end
    function f = get.available(MD)
        if isempty(MD.available), MD = refresh(MD); end; f = MD.available;
    end
    
    function source = get.source(MD), source = MD.data.Properties.CustomProperties.source; end
    function MD = set.source(MD,source)
        if ~isprop(MD.data.Properties.CustomProperties,'source')
            MD.data = addprop(MD.data,'source','variable');
        end
        MD.data.Properties.CustomProperties.source = source; 
    end
    function src = getsourceof(MD,fld), src = [MD.source{fieldidx(MD,fld)}]; end
    function MD = setsourceof(MD,fld,val)
    % Equivalent to MD.source{key} = val, which doesn't work inside the class because overloaded 
    % subsref is not called.

        fld = cellstr(fld);
        ib = fieldidx(MD,fld);
        try
            if isscalar(ib)
                val = {cellstr(val)};
            else
                validateattributes(val,{'cell'},{'numel',numel(ib)});
                for j = 1:numel(ib)
                    val{j} = cellstr(val{j});
                    assert(numel(val{j}) == size(MD.data.(fld{j})));
                end
            end
        catch
            error('Expecting VAL to be a cell of n-cellstrings, for each n-column FLD');
        end
        MD.source(ib) = val;
    end
    
    function ib = fieldidx(MD,fld)
        [ia,ib] = ismember(fld,MD.data.Properties.VariableNames);
        if ~all(ia)
            fld = cellstr(fld);
            error(shortliststr(fld(~ia),'Unknown field','quotes',''''));
        end
    end
    
    function MD = refresh(MD)
    % Update transient properties MD.dark & MD.available
    
        opt = MD.options;

        dt = MD.timestep;
        if isempty(dt) || isnan(dt), MD = checktimestamps(MD,false); end
        
        f = ~any(isfinite(MD.data{:,:}),2);
        if opt.fillgaps > 0
        % Don't flag small gaps as missing data
            ixm = find(~f);
            d = interp1(ixm,ixm,1:N,'next','extrap')-interp1(ixm,ixm,1:N,'previous','extrap') - 1; 
            f(d <= opt.fillgaps) = false;
        end
        MD.missing = f;
        
        FLD = intersect(fieldnames(MD),MeteoData.varnames('irradiance'));
        if ~isempty(FLD)
            dunkel = ~any(MD.data{:,FLD} > MD.options.minGHI,2);
        else
            dunkel = false(MD.Nt,1);
        end
        if isfield(MD,'sunel')
            dunkel = dunkel | MD.data.sunel < MD.options.minSunEl;
        end
        MD.dark = dunkel;
        
        f = intersect(MD.data.Properties.VariableNames,MeteoData.varnames('irradiance'));
        MD.available = any(isfinite(MD.data{:,f}),2);
    end
 
    function MD = set.interval(MD,k)
        if isempty(MD.interval), MD.interval = k; end
        [K,offset] = checksummarization({MD.interval,k});
        if ~any(K == 'i') && diff(offset) ~= 0
            MD.t = MD.t + diff(offset)*MD.timestep; %#ok<MCSUP>
        end
        MD.interval = K(2);
    end

    function x = isfield(MD,fld)
        if ~iscell(fld), fld = {fld}; end
        x = ismember(fld,fieldnames(MD));
        x(x) = ~cellfun(@isempty,fld(x));
    end
    
    function MD = rmfield(MD,fld)
        ib = fieldidx(MD,fld);
        MD.data(:,ib) = [];
        % MD.data.Properties.CustomProperties.source(ib) = [];
        MD.flags.data(:,ib) = [];
    end
    
    function [MD,X,F] = rmsource(MD,src)
    % [MD,X,F] = RMSOURCE(MD,SRC) remove data and flags for source(s) SRC. Removed data X and
    %   flags F are returned as second and third arguments.
    
        [X,F,ib,ic] = getbysource(MD,src);
        for j = unique(ib)
            inj = (ib == j);
            fld = MD.data.Properties.VariableNames{j};
            if isempty(setdiff(1:size(MD.data.(fld),2),ic(inj)))
                MD.data.(fld) = [];
                MD.flags.data.(fld) = [];
            else
                MD.data.(fld)(:,ic(inj)) = [];
                MD.data.Properties.CustomProperties.source{j}(ic(inj)) = [];
                MD.flags.data.(fld)(:,ic(inj)) = [];
            end
        end
    end
    
    function [X,F,var_ix,col_ix] = getbysource(MD,src)
    % [X,F,iv,ic] = GETBYSOURCE(MD,SRC) extract data X and flag F column(s) with source ID SRC
    %   return also variable- and column-indices iv, ic. 
    
        src = cellstr(src);
        if size(src,1) > 1, src = src'; end
        validateattributes(src,{'cell'},{'nonempty','vector'},'','src');
        try 
            [var_ix,col_ix] = cellfun(@(x) ismember(src,x),MD.source,'unif',0);
            var_ix = cat(1,var_ix{:}); col_ix = cat(1,col_ix{:});
            assert(all(any(var_ix,1)));
        catch
            MD = checksources(MD);
            [var_ix,col_ix] = cellfun(@(x) ismember(src,x),MD.source,'unif',0);
            var_ix = cat(1,var_ix{:}); col_ix = cat(1,col_ix{:});
            notfound = ~any(var_ix,1);
            assert(~any(notfound),'Failed to find %s',shortliststr(src(notfound),'source'));
        end
        col_ix = col_ix(var_ix);
        [var_ix,~] = find(var_ix);

        X = cell(1,numel(src));
        F = cell(1,numel(src));
        for j = 1:numel(src)
            fld = MD.data.Properties.VariableNames{var_ix(j)};
            X{j} = MD.data.(fld)(:,col_ix(j));
            F{j} = MD.flags.data.(fld)(:,col_ix(j));
        end
        X = cat(2,X{:});
        F = cat(2,F{:});
    end
    
    function MD = addfield(MD,fld,val,src,overwrite)
        Nc = size(val,2);
        val = compatiblesize(val,MD.t);
        if nargin < 4 || isempty(src), src = []; end
        if nargin < 5 || isempty(overwrite), overwrite = false; end
        
        if ischar(src) || isempty(src), src = {src}; end
        if Nc > 1 && isscalar(src), src = repmat(src,1,Nc); end
        validateattributes(src,'cell',{'vector','numel',Nc});
        
        if overwrite && isfield(MD,fld), MD = rmfield(MD,fld); end
        for j = 1:Nc
            MD = addsource(MD,fld,val(:,j),src{j});
        end
    end
    
    function [MD,ib,n] = addsource(MD,fld,val,src,overwrite)
    % MD = ADDSOURCE(MD,FLD,VAL,[SRC]) - Add column with source-ID SRC to existing field FLD,
    %   creating a new field if it doesn't exist. Just specifying FLD uses the default source
    %   ID 'FLD.N' where N is the number of columns of MD.FLD (after adding the new VAL).
    
        if nargin < 5 || isempty(overwrite), overwrite = false; end
        
        validateattributes(fld,{'char'},{'nonempty'},'','fld');
        % validateattributes(val,{'numeric','logical'},{'2d','size',[MD.Nt,NaN]},'','val');
        if nargin < 4 || isempty(src), src = []; 
        else
            if iscellstr(src) || isstring(src), src = char(src); end
            validateattributes(src,{'char'},{'nonempty'},'','src');
        end
        
        [~,ib] = ismember(fld,MD.data.Properties.VariableNames);

        if ib == 0
            n = 1;
            MD.data.(fld) = val;
            MD.flags.data.(fld) = zeros(MD.Nt,1,meteoQC.TYPE);
            ib = size(MD.data,2);
            overwrite = false;
        else
            n = size(MD.data.(fld),2)+1;
            if overwrite
                [overwrite,ic] = ismember(src,cellstr(MD.source{ib}));
                if overwrite, n = ic; end
            end
            MD.data.(fld)(:,n) = val;
            MD.flags.data.(fld)(:,n) = 0;
        end
        
        if overwrite, return; end
        
        if isempty(src), src = sprintf('%s.%d',fld,n); end

        srclist = [MD.data.Properties.CustomProperties.source{:}];
        usrc = matlab.lang.makeUniqueStrings(src,srclist);
        if ~strcmp(src,usrc)
           warning('Using %s instead of existing ID %s',usrc,src);
           src = usrc;
        end
        MD.data.Properties.CustomProperties.source{ib}(n) = {src};
    end
    
    function MD = fieldop(MD,fld,op)
        src = getsourceof(MD,fld);
        [X,F] = getbysource(MD,src);

        X = op(X);
        G = zeros(size(X),meteoQC.TYPE); 
        for k = 2:numel(MD.flags.flags)
            G = bitset(G,k,op(bitget(F,k)));
        end
        MD.data.(fld) = X;
        MD.flags.data.(fld) = G;
    end
    
    function MD = qcfield(MD,varargin)
    % MD = QCFIELD(MD,OKFCN,BIT/KEY,FLD,[FLAGIFNAN]) - shortcut for METEOQC.CHECKFIELD
        MD.flags = checkfield(MD.flags,MD,varargin{:});
    end
    
    function [MD,X,F,var_ix,col_ix] = qcsource(MD,okfcn,bit,src,flagifnan)
    % [MD,X,F,var_ix,col_ix] = QCSOURCE(MD,OKFCN,BIT/KEY,SRC,[FLAGIFNAN]) - shortcut for
    %   METEOQC.CHECKFIELD, but limited to a single variable (identified by source)
    
        narginchk(4,5);
        if nargin < 5, flagifnan = false; end
        
        [X,F,var_ix,col_ix] = getbysource(MD,src);
        fld = MD.data.Properties.VariableNames(var_ix);
        
        if ~isnumeric(bit), [bit,MD.flags] = flagbit(MD.flags,bit); end
        validateattributes(bit,{'numeric'},{'integer','scalar','positive'});
        validateattributes(okfcn,{'function_handle'},{'scalar'});
        
        for j = 1:numel(fld)
           F(:,j) = meteoQC.check(X(:,j),okfcn,bit,F(:,j),flagifnan);
           MD.flags.data.(fld{j})(:,col_ix) = F(:,j);
        end
    end
    
    MD = removemissing(MD,fld,E,W,quiet)
    
    function X = renamevars(X,varargin)
        X.data = renamevars(X.data,varargin{:});
        X.flags.data = renamevars(X.flags.data,varargin{:}); 
    end
    
    function MD = MeteoData(X,varargin)
    % obj = MeteoData(S,..) - parse table, timetable, or structure S into MeteoData object.
    %
    % obj = MeteoData(..,'prop',val) - set/override property (e.g. location, sensors, t,...)
    %
    % obj = MeteoData(..,'var','alias') - for variable names recognized by FINDFIELDS, tells the
    %       constructor to place field 'alias' of S as a column of obj.var, with source 'alias'.
    %       The same can be achieved with ..,'sensors',S, if one of the MeteoSensor objects in S
    %       has ID == 'alias' and type == 'var'.
    %
    % obj = MeteoData(..,'opt',val) - Set one of the options used by COMPLETEMETEODATA, LOADOBJ,
    %       REFRESH, etc.
    %
    % See also: METEODATA, COMPLETEMETEODATA, METEOQC, METEOSENSOR

    % meteo = 'clearskycheck', 'downsample', 'fillgaps', 'filter', 'interval', 'minCSfraction', 'separationmodel', 'tolerance', 'useweb'
    % outliers = 'P', 'beyondtol', 'error', 'warning', 'withintol';

        if nargin > 0
            validateattributes(X,{'MeteoData','table','timetable','struct'},{},'X');
            
            [opt,ADD,varargin] = parseargs(varargin{:});

            % See checksources for the use of Idx
            if ~isempty(ADD.source), Idx = {ADD.source}; else, Idx = {}; end
            
            if ~isempty(ADD.sensors)
            % Use sensors.ID and .type as dictionary to rename variables
                validateattributes(ADD.sensors,{'MeteoSensor'},{},'','sensors');
                [X,Idx{end+1}] = renamefields(X,[{ADD.sensors.ID}',{ADD.sensors.type}'],...
                    '-ignorecase','cat',2);
            end

            % Look for common variable aliases, 
            [X,Idx{end+1}] = MeteoData.findfields(X,varargin{:});
            Idx{end} = rmfield(Idx{end},intersect(fieldnames(Idx{end}),properties(MeteoData())));
            
            X = {X};
        else
            X = {};
        end
        
        % superclass constructor handles conversion from table/timetable/structure
        MD@table_ish(X{:});

        if nargin == 0, return; end
        
        if ~isempty(MD.options)
            [~,~,msg] = comparestruct(MD.options,opt,'and');
            if ~isempty(msg)
                warning('Overriding options:\n%s',strjoin(msg,newline())); 
            end
            [~,~,msg] = comparestruct(MD.options,opt,'diff');
            if ~isempty(msg)
                warning('Removing obsolete/unknown options:\n%s',strjoin(msg,newline())); 
            end
        end
        MD.options = opt;
        
        if ~isempty(MD.data.Properties.UserData) && isstruct(MD.data.Properties.UserData)
            R = MeteoData.findfields(MD.data.Properties.UserData);
            ADD = completestruct(ADD,R); % ,'valid',@(x) ~isempty(x));
        end

        varnames = MD.data.Properties.VariableNames;
        known = ismember(varnames,MeteoData.varnames('all'));
        if isempty(MD.data.Properties.VariableUnits)
           MD.data.Properties.VariableUnits = repmat({''},1,size(MD.data,2)); 
        end
        if isempty(MD.data.Properties.VariableDescriptions)
           MD.data.Properties.VariableDescriptions = repmat({''},1,size(MD.data,2)); 
        end
        for j = 1:numel(known)
            if ~known(j), continue; end
            [d,u] = MeteoData.whatis(varnames{j});
            if isempty(MD.data.Properties.VariableDescriptions{j}) || ...
                    strcmp(MD.data.Properties.VariableDescriptions{j},varnames{j})
                MD.data.Properties.VariableDescriptions{j} = d;
            end
            if isempty(MD.data.Properties.VariableUnits{j})
                MD.data.Properties.VariableUnits{j} = u;
            end
        end
        
        % Idx might, at this point, hold different "layers" of renaming steps,
        % Merge them into a single map: fieldnames(MD) --> MD.source (IDs)
        MD = checksources(MD,Idx{:});

        if ~isprop(MD.data.Properties.CustomProperties,'sensors')
           MD.data = addprop(MD.data,'sensors','variable'); 
        end
        % [MD,notmatched] = getsensordata(MD,sensors)

        Loc = parselocation(ADD,'optional','all');
        locfields = {'latitude','longitude','altitude','TimeZone','name'};
        Loc = rmfield(Loc,setdiff(fieldnames(Loc),locfields));
        ADD = rmfield(ADD,fieldnames(Loc));
        
        if isempty(ADD.location), ADD.location = struct(); end
        ADD.location = completestruct(ADD.location,Loc);
        ADD.location = parselocation(ADD.location,'optional','all');
            
        for f = {'location','options'}
            fld = f{1};
            if isempty(ADD.(fld)), continue; end
            if isempty(MD.(fld)), MD.(fld) = struct(); end
            [~,~,msg] = comparestruct(ADD.(fld),MD.(fld),'and');
            if ~isempty(msg)
                warning('Overriding %s fields:\n%s',fld,msg);
            end
            MD.(fld) = completestruct(ADD.(fld),MD.(fld));
        end
        if ~isempty(MD.location)
        % Parse location: guess TimeZon based on longitude, and set altitude = 0, if necessary
            MD.location = parselocation(MD.location,'useweb',opt.useweb,'-soft');
        end
        
        % Initialize (and merge) MeteoQC flags object
        MD.flags = meteoQC(MD,ADD.flags);
        
        % Parse sensor details, if available (don't assume defaults, yet)
        MD = checksensors(MD,ADD.sensors); 
        
        if ~isempty(ADD.info)
            ADD.info = cellstr(ADD.info);
            MD.info = uniquecell(cat(1,MD.info(:),ADD.info(:)));
        end
             
        for f = {'t','timestep','interval'}
            fld = f{1};
            if isempty(ADD.(fld)), continue; end
            if isempty(MD.(fld)) || ~any(isfinite(MD.(fld))), MD.(fld) = ADD.(fld); continue; end
            if ~isequal(MD.(fld),ADD.(fld))
                error('Incompatible (existing) MD.%s and provided value',fld);
            end
        end
        
        used = {'t','timestep','interval','sensors','flags','info',...
                'location','options','source','timestep','interval'};
            
        % still to handle: 'models','status',
            
        ADD = rmfield(ADD,intersect(fieldnames(ADD),used));
        MD.data.Properties.UserData = ADD;
                
        % OPT = {'sensors','info','uncertainty'};

        % Ensure unique, valid timestamps
        MD = checktimestamps(MD);
    end
    
    MD = checksensors(MD,sensors,assumedefaults)

    function MD = checktimestamps(MD,regular,dt)
    % obj.checktimestamps() - Timestep parsing & regularization (assume UTC, if missing)

        if nargin < 2, regular = true; end
        if nargin < 3, dt = []; end

        if isempty(MD.interval)
            try
                [u,dt,idx,ia] = parsetime(MD.t,'step',dt,'interval',MD.interval,'-regular');
                warning('Assuming center-of-interval summarization');
                MD.interval = 'c';          
            catch ERR
                if ~isempty(dt) || ~strcmp(ERR.identifier,'parsetime:notgridded')
                    rethrow(ERR); 
                end
                [u,~,idx,ia] = parsetime(MD.t,'interval','i','-sorted','-unique');
                MD.interval = 'i';
                dt = NaN;
            end
        else
            MD.interval = checksummarization(MD.interval);
            if MD.interval(end) == 'i'
                [u,~,idx,ia] = parsetime(MD.t,'step',dt,'interval',MD.interval,'-sorted','-unique');
                dt = NaN;
            else
                [u,dt,idx,ia] = parsetime(MD.t,'step',dt,'interval',MD.interval,'-regular');
            end
        end

        MD.t = u(idx(ia));
        MD.timestep = dt;
        MD.interval = MD.interval(end);

        notthere = true(numel(u),1);
        notthere(idx) = false;

        n = accumarray(ia,1);
        % repeated = idx(n > 1);
            
        if any(n > 1) || numel(u) > MD.Nt || ~issorted(idx(ia))
        % Average repeated time-steps, re-sort and regularize

            F = sparse(idx(ia),1:numel(ia),1./n(ia),numel(u),numel(ia));
            % F = F + sparse(find(notthere),1,NaN,numel(u),numel(ia));
            % MD = filterstructure(MD,F);
            
            MD.data = filterstructure(MD.data,F);
            MD.flags = filterstructure(MD.flags,F);
            if ~isempty(MD.uncertainty)
                MD.uncertainty.data = filterstructure(MD.uncertainty.data,F);
            end
            [MD.t,dt] = parsetime(MD.t,'-gridded');
            if (MD.interval ~= 'i'), MD.timestep = dt; end
        
            if regular
                MD.data{notthere,:} = NaN(size(MD.data{notthere,:}));
                MD.t = u;
            else
                MD.data(notthere,:) = [];
                MD.flags.data(notthere,:) = [];
                if ~isempty(MD.uncertainty)
                    MD.uncertainty.data(notthere,:) = [];
                end
            end
            MD = refresh(MD);
        end
    end
    
    function MD = checkUTCoffset(MD,varargin)
    % obj.checkUTCoffset(...) - UTC-offset check (and correction), wrapper for CHECKUTCOFFSET.
    %   Any arguments are passed as flags / name-value pairs to CHECKUTCOFFSET.
    
        parselocation(MD.location);
        if ~isregular(MD.data), MD = checktimestamps(MD); end
        try
            if isfield(MD,'GHI')
                GHI = mean(MD.data.GHI,2,'omitnan');
                offset = checkUTCoffset(GHI,MD.t,MD.location,'interval',MD.interval,varargin{:}); 
            elseif isfield(MD,'BNI')
                BNI = mean(MD.data.GHI,2,'omitnan');
                MD.interval = 'c';
                [~,CSBNI] = pvlmod_clearsky_ineichen(MD.location,MD.t);
                offset = checkUTCoffset(BNI,MD.t,CSBNI,varargin{:}); 
            else
                error('Cannot check UTC offset without GHI/BNI and location/clear-sky');
            end
        catch ERR
            warning('meteodata:checkutc','Failed to check UTC offset: %s',ERR.message);
            return;
        end        
        if offset~=0
            MD.t = MD.t - offset;
            MD = getsunpos(MD,0,1);
            MD = refresh(MD);
        end
    end

    % [MD,B] = bestimate(MD,varargin)

    MD = checksources(MD,varargin)
    
    MD = filter(MD,varargin)
    
    function MD = filterstructure(MD,varargin)
    % MD = FILTERSTRUCTURE(MD,...) - overloaded FILTERSTRUCTURE method for filtering, resampling,
    %   or downsampling of METEODATA objects.
        
        original.interval = MD.interval;
        original.timestep = MD.timestep;
        original.TimeZone = MD.t.TimeZone;
        
        if ~any(MD.interval == 'ic'), MD.interval = 'c'; end
        MD.t.TimeZone = 'UTC';

        MD.data = filterstructure(MD.data,varargin{:});
        MD.flags = filterstructure(MD.flags,varargin{:});
        if ~isempty(MD.uncertainty)
            MD.uncertainty.data = filterstructure(MD.uncertainty.data,varargin{:});
        end
        MD = checktimestamps(MD);
        
        if ~isequal(MD.timestep,original.timestep)
        % Update effective solar position
        % TODO: let solar position be irradiance weighted, i.e. try to keep closure
        
            MD = refresh(MD);
            fld = {'sunel','sunaz','hourangle'};
            fld = fld(isfield(MD,fld));
            if ~isempty(fld)
                MD = rmfield(MD,fld);
                MD = getsunpos(MD,false);
            end
        end

        MD.interval = original.interval;
        MD.t.TimeZone = original.TimeZone;
        MD = refresh(MD);
    end

    function C = horzcat(A,varargin)
    % C = HORZCAT(A,B,..), C = [A,B,..] - horizontal concatenation (add columns for METEODATA 
    %   objects with common timestamps and the same location.
        
        C = A;
        switch nargin
        case 1, return;
        case 2, B = varargin{1}; % continue below
        otherwise
        % PROVISIONAL: code below should be rewritten to take several arguments  
            for j = 1:numel(varargin)
                C = horzcat(C,varargin{j});
            end
            return;
        end
        
        EQ = {'interval','timestep','Nt','t','location'};
        for j = 1:numel(EQ)
            if isempty(A.(EQ{j})), C.(EQ{j}) = B.(EQ{j});
            elseif ~isempty(B.(EQ{j})) && ~isequal(A.(EQ{j}),B.(EQ{j}))
                error('Expecting objects with common %s',EQ{j});
            end
        end

        allflags = unique([A.flags.flags,B.flags.flags],'stable');
        C.flags = setflags(A.flags,allflags);
        B.flags = setflags(B.flags,allflags);
        
        newfld = setdiff(fieldnames(B),fieldnames(A));
        if ~isempty(newfld)
            C.data = [C.data,B.data(:,newfld)];
            C.flags.data = [C.flags.data,B.flags.data(:,newfld)];
        end
        
        [ib,ia] = ismember(fieldnames(B),fieldnames(A));
        ia = ia(ib);
        ib = find(ib);
        for j = 1:numel(ib)
            fld = A.data.Properties.VariableNames{ia(j)};
            C.data.(fld) = [A.data{:,ia(j)},B.data{:,ib(j)}];
            C.flags.data.(fld) = [A.flags.data{:,ia(j)},B.flags.data{:,ib(j)}];
            C.source{ia(j)} = [A.source{ia(j)},B.source{ib(j)}];
        end
        
        C.sensors = [A.sensors,B.sensors];
        C.models = [A.models,B.models];
    end
    
    function C = vertcat(varargin)
    % C = VERTCAT(A,B,..), C = [A;B;..] - vertical concatenation (add time-steps for METEODATA 
    %   objects with common fields and sources.
        
        C = varargin{1};
        if nargin == 1, return; end
        
        % x = cat(1,varargin{:}); 
        
        EQ = {'interval','timestep','location'};
        for j = 1:numel(EQ)
            for k = 1:numel(varargin)
                B = varargin{k};
                if isempty(C.(EQ{j})), C.(EQ{j}) = B.(EQ{j});
                elseif ~isempty(B.(EQ{j})) && ~isequal(C.(EQ{j}),B.(EQ{j}))
                    error('Expecting objects with common %s',EQ{j});
                end
            end
        end

        allflags = cellfun(@(x) x.flags.flags,varargin,'unif',0);
        allflags = unique([allflags{:}],'stable');
                
        allsources = cellfun(@(x) [x.source{:}],varargin,'unif',0);
        [allsources,ic] = unique([allsources{:}],'stable');
        
        allfields = cellfun(@(x) repelem(x.fieldnames,cellfun(@numel,x.source)),varargin,'unif',0);
        allfields = cat(1,allfields{:});
        allfields = allfields(ic);
        
        for k = 1:numel(varargin)
            A = varargin{k};
            A.flags = setflags(A.flags,allflags);
            ic = ~ismember(allsources,[A.source{:}]);
            if any(ic)
                for j = find(ic)
                    A = A.addfield(allfields{j}, NaN, allsources{j}, false);
                end
            end
            if k == 1
                C = A;
                fld = C.fieldnames;
            else
                [ia,ib] = ismember(fld,A.fieldnames);
                assert(all(ia),'Something went wrong');
                if ~isequal(ib',1:numel(ia))
                    A.data = A.data(:,ib);
                    A.flags.data = A.flags.data(:,ib);
                end
                for j = 1:numel(fld)
                    [ia,ib] = ismember(C.source{j},A.source{j});
                    assert(all(ia),'Something went wrong');
                    if ~isequal(ib,1:numel(ia))
                        A.data.(fld{j}) = A.data.(fld{j})(:,ib);
                        A.flags.data.(fld{j}) = A.flags.data.(fld{j})(:,ib);
                    end
                end
            end
            varargin{k} = A;
        end
        TT = cellfun(@(x) x.data,varargin,'unif',0);
        C.data = cat(1,TT{:});
        TT = cellfun(@(x) x.flags.data,varargin,'unif',0);
        C.flags.data = cat(1,TT{:});
        
        allsensors = cellfun(@(x) x.sensors,varargin,'unif',0);
        allsensors = cat(1,allsensors{:})';
        [~,ic] = uniquecell({allsensors.ID},'stable');
        C.sensors = allsensors(ic);
        
        C = checktimestamps(C,false,C.timestep);
    end
    function x = cat(d,varargin) 
        validateattributes(d,{'numeric'},{'integer','scalar','positive','<',3});
        if d == 1, x = vertcat(varargin{:}); else, x = horzcat(varargin{:}); end
    end
          
    [UNC,SENS] = getuncertainty(MD,tolerance)
    
    varargout = plot(MD,type,varargin)
    
    function MD = cleanup(MD)
    % MD = CLEANUP(MD) - ensure that ancilliary variables [ MD.varnames({'ambient','clearsky'}) ]
    %   are vector fields. Remove columns copied from MERRA2 data (unless no other channel is
    %   available) and pick a single clear-sky model(FITCLEARSKY > MCCLEAR > INEICHEN).
        
        MD = checksources(MD);
        allfields = MD.data.Properties.VariableNames;
        
        % remove reduntant merra2 fields (already used for gap-filling)
        fld = intersect(allfields,MeteoData.varnames('ambient'));
        for k = 1:numel(fld)
            j = fieldidx(MD,fld{k});
            
            if numel(MD.source{j}) <= 1, continue; end
            
            idx = contains(MD.source{j},'merra2');
            if any(idx) && ~all(idx)
                MD = rmsource(MD,MD.source{j}(idx));
            end
            
            if numel(MD.source{j}) > 1
                warning('cleanup:average','Using simple average of %s',shortliststr(MD.source{j}));
                % MD.data.(fld{k}) = mean(MD.data{:,j},2,'omitnan');
                MD = MD.fieldop(fld{k},@(X) mean(X,2,'omitnan'));
            end
        end
        
        PRIORITY = {'fitclearsky','mcclear','linketurbidity'};
        
        % remove reduntant clear-sky fields (already used for gap-filling)
        fld = intersect(allfields,MeteoData.varnames('clearsky'));
        for k = 1:numel(fld)
            j = fieldidx(MD,fld{k});
            
            if numel(MD.source{j}) <= 1, continue; end

            for s = 1:numel(PRIORITY)
                idx = contains(MD.source{j},PRIORITY{s});
                if any(idx) && ~all(idx)
                    MD = rmsource(MD,MD.source{j}(~idx));
                end
                if numel(MD.source{j}) <= 1, break; end
            end
            
            if numel(MD.source{j}) > 1
                warning('cleanup:average','Using simple average of %s',shortliststr(MD.source{j}));
                % MD.data.(fld{k}) = mean(MD.data{:,j},2,'omitnan');
                MD = MD.fieldop(fld{k},@(X) mean(X,2,'omitnan'));
            end
        end
    end
    
    function [B,P,U] = legacy(MD)
    % [B,P,U] = LEGACY(MD) - generate simplified structures, compatible with older code versions:
    %
    %   B: structure with vector fields (single-channel) GHI, DHI, BNI, Ta, vw...
    %   P: structure with vector fields (solar position variables)
    %   U: irradiance uncertainty
    
        [MD,B,U] = bestimate(MD);
        % B = timetable2struct(MD.data(:,{'kn','kt','kd'}));
        MD = cleanup(MD);

        allfields = MD.data.Properties.VariableNames;
        
        f = intersect(allfields,MeteoData.varnames('sunpos'));
        P = table2struct(MD.data(:,f),'toscalar',true);
        P = renamefields(P,{'sunel','El';'sunaz','Az';'hourangle','w';'declination','dec'});

        f = MeteoData.varnames({'ambient','clearsky','indices'});
        f = intersect(allfields,f);
        B = completestruct(B,table2struct(MD.data(:,f),'toscalar',true));
        
        B.ENI = P.ENI;
        P = rmfield(P,'ENI');
    end
    
end
methods (Static = true)
%     function MD = loadobj(MD)
%        MD = refresh(MD); % does not work :( since MD is not an object yet
%     end
    
    function MD = import(filename,sensorsfile,filetype,varargin)
    % MD = METEODATA.IMPORT(filename,sensorsfile,filetype,...) - wrapper for METEODATA.GETMETEODATA
    %   to include a configuration (*.sensors) file. 
    
        if nargin < 1 || isempty(filename), filename = pickfile({'*.meteo'}); end
        if nargin < 2 || (isempty(sensorsfile) && ~iscell(sensorsfile))
            sensorsfile = pickfile('*.sensors',Inf,'ui',1); 
        else
            sensorsfile = cellstr(sensorsfile);
        end
        if nargin < 3 || isempty(filetype), filetype = ''; end
        
        sens = cellfun(@MeteoSensor.readsensorspecs,sensorsfile,'unif',0);
        sens = cat(1,sens{:});
        
        [S, time, Loc] = MeteoData.getMeteoData(filename,filetype);
        MD = MeteoData(S,'t',time,'location',Loc,'sensors',sens,varargin{:});
    end
    
    [S, time, Loc] = getMeteoData(filename,filetype);
    
    function f = varnames(types)
    % C = METEODATA.VARNAMES(TYPES) - Return named subsets of known MeteoData properties, e.g.
    %
    %   MeteoData.varnames('irradiance') = {'DHI', 'BNI', 'GHI', 'GTI', 'USW'}
    %   MeteoData.varnames('all') == fieldnames(MeteoData())
    %
    %   setxor(fieldnames(MeteoData),MeteoData.varnames('all')) should be empty!
    
        if nargin == 0, types = 'numeric'; end
        
        p = {'meta',{'location', 'info', 'options', 'sensors', 'uncertainty', 'flags', 'var'};
            'time',{'t', 'timestep', 'interval'};
            'irradiance',{'DHI', 'BNI', 'GHI', 'GTI', 'USW'};
            'sunpos',{'ENI', 'sunel', 'sunaz', 'hourangle', 'declination'};
            'ambient',{'Patm', 'RH', 'tpw', 'windir', 'vw', 'Ta', 'soiling'};
            'clearsky',{'clearsky', 'CSGHI', 'CSBNI', 'CSDHI', 'TL', 'AMa', 'AOD'};
            'indices',{'kd', 'kt', 'kn', 'albedo'};
            'dependent',{'Nt', 'dark', 'available'}};

        redundant = size(p,1);
        % p(end+1,:) = {'all',cat(2,p{1:end})};
        p(end+1,:) = {'numeric',cat(2,p{3:7,2})}; % irradiance,..,indices
        p(end+1,:) = {'measured',cat(2,p{[3,5],2})}; % irradiance & ambient

        [~,idx] = parselist(types,p(:,1),'key');
        f = cat(2,p{idx,2})';
        if any(idx > redundant), f = unique(f,'stable'); end
        
    end
    
    varargout = findfields(S,varargin)
    
    [s,u] = whatis(x,varargin)
    
    [xc,varargout] = fillgaps(X,dt,timescale,varname,excluded,varargin)
    
end

end
    
function [opt,prop,rest] = parseargs(varargin)
% Parse possible list of options destined for different purposes:
%   OPT -> MD.options, later used by loadobj and meteoQC.test
%   PROP direct MD.x property attributions (used by constructor)
%   REST - possibly a list of alias to be parsed by FINDFIELDS

    FLAGS = {'-verbose','-stickler','-useweb'};
    PARAM = {'verbose','useweb','stickler','minSunEl','minGHI','outliers'}; % + meteo.* fields
    PROPS = {'flags','info','location','models','options','interval',...
                'sensors','source','status','t','timestep'};
            
    opt = getSimOption('meteo');
    opt = completestruct(opt,getSimOption(PARAM));
    [opt,names] = nestedstruct2cell(opt);
    opt = cell2struct(opt,strrep(names,'outliers.','outlier_'));
    
    [opt,varargin] = parseoptions(varargin,FLAGS,opt);

    % Parse options
    parsestruct(opt,{'outlier_P','outlier_warning','outlier_error',...
                     'minCSfraction','minGHI','minSunEl','fillgaps'},'-r','-p','-s');
    parsestruct(opt,{'verbose','useweb','stickler'},'-s','-l');

    [prop,varargin] = getpairedoptions(varargin,PROPS,{[]});

    % % Two kinds of name-value pairs:
    % % 'location','Loc', ..  as in: search for 'location' in field 'Loc' (for FINDFIELDS)
    % % 'location', Loc, ..   as in: 'location' is object Loc
    % weresources = cellfun(@(f) ischar(prop.(f)),PROPS);
    % if any(weresources)
    %     prop = struct2cell(prop);
    %     args = [PROPS(weresources);prop(weresources)]';
    %     varargin = [varargin,args(:)'];
    %     prop(weresources) = {[]};
    %     prop = cell2struct(prop,PROPS);
    % end

    [opt.QC,rest] = parseoptions(varargin,{'-independent'},getSimOption('meteoQC'));
end
