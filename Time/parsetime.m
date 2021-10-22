function [T,dt,idx,ia] = parsetime(X,varargin)
% T = PARSETIME(X) - Convert whatever is in X {DATETIME, DATEVEC, PVL_MAKETIMESTRUCT, DATESTR,
%   DATENUM*...} into a DATETIME vector. Unless Time-Zone/UTCoffset information is included,
%   result will have T.TimeZone = 'UTC' (see 'TimeZone' option below for more details).
%
%   (*) Other numeric types [POSIX-time (s), POSIX-time (ms) and Excel] might be recognized based
%   on the scale of the input array, and converted with a warning.
%
% T = PARSETIME(..,'output',KEY) for KEY in {'datetime' (default),'datevec','datenum','struct'} 
%   selects a different output format. Numeric formats will be be adjusted to requested 'TimeZone'
%   only if provided explicitly (see below). 
%
%   NOTE: For compatibility with PVLIB. 'struct' is equivalent to PVL_MAKETIMESTRUCT.
%
% T = PARSETIME(..,'InputFormat',FMT) specify format for conversion of cellstr and charr-arrays. 
% T = PARSETIME(..,'ConvertFrom',S,..) just calls DATETIME with the same arguments
%
% T = PARSETIME(..,'TimeZone',TZ) specify Time-Zone of output T. TZ can be anything recognized
%     by DATETIME (e.g. 'local','UTC',integers,'+08:00','Europe/Berlin',...) or the keywords
%     {'keep','none'}. This can have diverse effects:
%
%   a) if X already contains Time-Zone information (e.g. DATETIME, or PVL_MAKETIMESTRUCT), the
%     resulting T will be X converted to the new Time-Zone (i.e. time offset might be applied).
%
%   b) if X has no Time-Zone information, the resulting T quietly assumes new TZ, except if TZ
%     is empty (default!), in this case a warning will be issued, and 'UTC' will be assumed. 
%
%   c) 'TimeZone','keep',.. overrides the behavior (b), T will keep whatever time-zone information  
%      is included in X, including T.TimeZone = '' (if that is the case).
%
%   d) 'TimeZone','none',.. will delete any time-zone information included in X.
%
% [T,dT] = PARSETIME(X,...) - attempt to detect time-step dT (minutes) automatically, setting
%   dT = mode(round(diff(T)/U))·U. Setting grid-units U to 1 minute, or to 1 second (if using
%   1 minute results in dT = 0), and ignoring repeated timesteps (diff(T) < U/2).
%
% [T,dT] = PARSETIME(..,'gridunits',U) - where U is a DURATION scalar, or a double [minutes],
%   forces dT to be a multiple of U, and matches the output T to the closest points in a regular 
%   grid with timestep U.
%
% [T,dT] = PARSETIME(..,'-gridded') - equivalent to using 'gridunits',dT. NOTE: '-gridded' is set 
%   as default whenever 'gridunits' ~= [], 'latticeunit' ~= [], or when nargout > 1. 
%   Use ..,'gridded',false to override this behaviour.
%
% [T,~] = PARSETIME(X,'step',DT) - for scalar duration DT, override search for time-step dT.
%   Settings for 'gridunits' will be ignored.
%
% [T,dT] = PARSETIME(..'step',@F) - use @F e.g. @MIN, @MEDIAN instead of MODE to find DT.
%   Consider the series X = datetime() + minutes([0,1,3,5,7]). PARSETIME(X,'-gridded') would
%   find a dT of 2 min, and force T(1) = X(1) - minutes(1), to fit into the grid. Using @min,
%   on the other hand, will find a 1-min grid with missing values on T(1) + 2,4 and 6 min.
%
% [S,dT,IDX] = PARSETIME(..,'-sorted') - return a sorted permutation of T, along with the inverse
%   sort index IDX, such that S(IDX) = T. Using '-sorted' with NARGOUT < 3 will throw an error 
%   if the points were not originally sorted.
%
% [R,dT,IDX] = PARSETIME(..,'-regular') - return a regular (sorted and equally spaced) time-step
%   grid R and an index vector IDX, such that all |R(IDX) - T| < |F·dT|. Note that, in the
%   general case, numel(R) >= numel(T), i.e. not all grid points must exist in T. However, using
%   '-regular' with NARGOUT < 3 will throw an error unless the points are sorted and without gaps.
%
% [U,dT,IDX] = PARSETIME(..,'-unique') - check for uniqueness, and issue an error if any time-
%   stamps are repeated [after fitting to grid]. 
%
% [U,dT,IC,IA] = PARSETIME(..,['-unique' implicit]) - Allow repeated timesteps, with a warning.
%   In this case, IDX will be reduced to unique position indices IC, and multiplicity will be 
%   included in IA, that is: IDX = IC(IA) for both '-sorted' and '-regular'.
%
% [...] = PARSETIME(..,'latticeunit',S) for DURATION scalar S, forces the first grid point to be an
%   integer fraction of S (since the start of the year). By default, S = seconds(1). Notice that
%   for units like S = days(1) or S = years(1), the result might change depending on 'TimeZone'.
%
% [...] = PARSETIME(..,'interval',K) for K in {'i','b','c','e'} implying instantaneous, beginning-
%   center- and end-of-interval (resp.), introduces an offset in the "anchor point" of the time-
%   steps, so that 'latticeunit' always matches with the beginning of the (finite) interval. e.g.
%   for dt = 10 minutes:
%      PARSETIME(...,'latticeunit',hour(1),'interval','b') will use 00:00, 00:10, 00:20 ...
%      PARSETIME(...,'latticeunit',hour(1),'interval','c') will use 00:05, 00:15, 00:25 ...
%      PARSETIME(...,'latticeunit',hour(1),'interval','e') will use 00:10, 00:20, 00:30 ...
%
% [...] = PARSETIME(..,'interval',{K,L}) - has the effect described above, but it also changes the
%   summarization convention of the resulting T to L. E.g. 'interval',{'b',e'} will effectively
%   return T = X + dT, whereas 'interval',{'c',b'} will return T = X - dT/2.
%
% [...] = PARSETIME(..,'tolerance',F) - (default 0.5) restrict adjustments to the grid to |F·U|, 
%   throwing an error if any |X - T| > |F·U|.
%
% NOTE: depending on the input/output format, conversion errors of order ~10·eps(now()) might be
%   expected.
%
% TODO: turn this into a wrapper class for DATETIME, with properties 'regular','interval', etc.
%   to avoid repeated parsing
%
% See also: RETIME, TIMETABLE

    if nargin == 0, test(); return; end  % DEBUG

    [OPT,varargin] = getflagoptions(varargin,{'-gridded','-regular','-sorted','-unique'});
    if ~OPT.gridded, OPT.gridded = []; end
    % if ~OPT.regular, OPT.regular = nargout > 2; end
    OPT.unique = OPT.unique || nargout > 3;
    OPT.output = 'datetime';
    OPT.gridunits = [];
    OPT.tolerance = 0.5;
    OPT.epoch = [];
    OPT.latticeunit = [];
    OPT.step = [];
    OPT.interval = [];
    OPT.TimeZone = '';
    OPT.InputFormat = '';
    OPT.convertfrom = '';
    OPT = getpairedoptions(varargin,OPT,'restchk');
    
    if OPT.regular, OPT.gridded = true; end
    if isempty(OPT.gridded)
        OPT.gridded = nargout > 1 || ~isempty(OPT.gridunits) || ~isempty(OPT.latticeunit) || ...
                      ~isempty(OPT.step); 
    end
    if isempty(OPT.latticeunit), OPT.latticeunit = seconds(1); end
    assert(isa(OPT.latticeunit,'duration'), 'Expecting DURATION latticeunit');
    if isa(OPT.gridunits,'duration'), OPT.gridunits = minutes(OPT.gridunits); end

    if ~isempty(OPT.convertfrom), OPT.InputFormat = OPT.convertfrom; end
    if isempty(OPT.interval), OPT.interval = 'i'; end
    
    % Bring everything to a DATETIME vector...
    t = any2datetime(X,OPT.InputFormat,'epoch',OPT.epoch);
    
    
    % Parse time zone, including keyworkds {'keep','none',''}
    [OPT.TimeZone,timezone_key] = checktimezone(OPT.TimeZone);

    if timezone_key ~= 'K'
        % Make sure we have a TimeZone by now, issue warning if using default (UTC)
        if isempty(t.TimeZone)
            if isempty(OPT.TimeZone)
                warning('parsetime:notimezone','Assuming UTC time-zone');
                t.TimeZone = 'UTC';
            else
                t.TimeZone = OPT.TimeZone;
            end
        else
            t.TimeZone = checktimezone(t.TimeZone,'UTC');
        end        
        
        % WORK INTERNALLY WITH UTC: DATETIME.UNIQUE does not detect DST jumps!
        if isempty(OPT.TimeZone), OPT.TimeZone = t.TimeZone; end
        if ~isequal(t.TimeZone,'UTC'), t.TimeZone = 'UTC'; end
    end

    sz = size(t);
    t = t(:);
    
    % Parse interval summarization, calculate offsets
    [OPT.interval,lattice_offset] = checksummarization(OPT.interval);
    switch numel(OPT.interval)
        case 1
            lattice_offset(2) = lattice_offset;
        case 2
            if OPT.interval(1) == 'i' && ~any(OPT.interval(2) == 'ci')
                warning('Ignoring summarization offset: instantaneous measurements');
                lattice_offset(:) = 0;
            elseif OPT.interval(2) == 'i' && OPT.interval(1) ~= 'i'
                warning('Requested ''instant'' summarization: returning center of intervals');
                lattice_offset(2) = 0.5;
            end
            if diff(lattice_offset) ~= 0, OPT.gridded = true; end % we need dT!
        otherwise, error('Expecting scalar or 2-cell interval summarization');
    end
    
    if isa(OPT.step,'duration')
        dt = OPT.step;
        knownstep = dt > 0;
        if ~knownstep
            OPT.step = @mode;
            warning('Ignoring invalid override-step'); 
        end
    else
        if isempty(OPT.step), OPT.step = @mode; end
        assert(isa(OPT.step,'function_handle'),'Expecting DURATION or function-handle');
        dt = minutes(NaN);
        knownstep = numel(t) <= 1;
    end
    
    if OPT.regular && nargout < 3,  OPT.sorted = true; end
    if OPT.sorted || OPT.gridded || ~knownstep
        if ~issorted(t)
            [t,sidx] = sort(t);
        else
            sidx = [];
        end
        if OPT.sorted && nargout < 3 && ~isempty(sidx)
           error('parsetime:notsorted','Time steps are not sorted'); 
        end
    else
        sidx = []; 
    end

    if OPT.gridded && ~knownstep
    % Estimate dt round to OPT.gridunits or to minutes or seconds
        
        delta = diff(t);
        guessdt = @(U) OPT.step(round(delta(delta > U/2)/U))*U;
        
        if isempty(OPT.gridunits)
            dt = guessdt(minutes(1));
            if ~(dt > 0)
                dt = guessdt(seconds(1));
            end
        else
            dt = guessdt(OPT.gridunits);
        end
        knownstep = dt > 0;
        if ~knownstep && ~OPT.gridded
            warning('parsetime:step','Invalid step');
        end
    end  
    
    if OPT.gridded
    % Return a uniform time-grid, and position index idx
    
        U = dt;
        assert(knownstep,'parsetime:step','Invalid step');
    
        t = t - lattice_offset(1)*dt;

        if isempty(t)
            idx = [];
        elseif isscalar(t)
            idx = 1;
            b = dateshift(t,'start','year');
            t = b + round((t - b)/OPT.latticeunit)*OPT.latticeunit;
        else
            % Find a starting lattice point
            offset = fminbnd(@(x) rssq(mod((t - t(1))/U + 0.5 + x,1)-0.5),-0.5,0.5);
            t0 = t(1) - minutes(offset);
            
            % Make sure starting point is a multiple of OPT.latticeunit
            b = dateshift(t0,'start','year');
            t0 = b + round((t0 - b)/OPT.latticeunit)*OPT.latticeunit;
            if t0 - t(1) > OPT.tolerance*U, t0 = t0 - max(U,OPT.latticeunit); end
            
            tE = t(end)+OPT.tolerance*U;
            if ~isempty(sidx), t(sidx)= t; end % revert sort

            for j = 1:3

                r = (t0:U:tE)';
                [~,idx] = ismembertol(datenum(t),datenum(r),OPT.tolerance*days(U),'DataScale',1);
                
                % fine tune starting point...
                ok = idx > 0;
                offset = mode(t(ok)-r(idx(ok)));
                if offset == 0, break; end

                last_t0 = t0;
                t0 = t0 + offset;
                t0 = b + round((t0 - b)/OPT.latticeunit)*OPT.latticeunit;
                if t0 == last_t0, break; end
            end
            assert(all(idx > 0),'parsetime:notgridded',...
                   'Time steps are not gridded within tolerance!');

            if ~OPT.regular
                t = r(idx);
                idx = (1:numel(t))';
                if OPT.sorted && ~isempty(sidx)
                    t = t(sidx);
                    idx(sidx) = idx;
                end
            else
                t = r;
                if nargout < 3
                    if isempty(sidx), sidx = (1:numel(t))'; end
                    assert(max(sidx) <= numel(idx) && isequal(idx(sidx)',1:numel(t)),...
                        'parsetime:notregular','Time steps are not regular within tolerance');
                end
            end
            
            t = t + lattice_offset(2)*dt;
        end
    else
        idx = (1:numel(t))';
        if ~isempty(sidx)
            if OPT.sorted, idx(sidx) = idx; else, t(sidx)= t; end % revert sort, or adjust index
        end
    end
    
    if OPT.unique
        if OPT.regular
            [idx,~,ia] = unique(idx,'stable');
        else
            [s,~,ia] = unique(t(idx),'stable');
            idx = (1:numel(s))';
            if OPT.sorted
                [t,is] = sort(s);
                idx(is) = idx;
            else
                t = s;
            end
        end

        repeated = accumarray(ia,1) > 1; % bool. index on ic
        if any(repeated)
            
            repeated = find(repeated); 
            n = numel(repeated);
            repeated(10:end) = [];
            msg = arrayfun(@(j) [datestr(t(idx(j))) ', idx: ' shortliststr(find(ia == j))],repeated,'unif',0);
            msg = sprintf('Found %d repeated %s\n',n,...
                    shortliststr(msg,'time-step',8,'newlines',true));
            if nargout > 3
                warning('cmd:reptimestamps',msg); %#ok<SPWRN>
            else
                error('cmd:reptimestamps',msg); %#ok<SPERR>
            end          
        end
    else
       ia = (1:numel(idx))'; 
    end
    
    % Bring back to requested TimeZone
    if timezone_key ~= 'K' && ~isequal(t.TimeZone,OPT.TimeZone)
        t.TimeZone = OPT.TimeZone;
    end
    
    switch lower(OPT.output)
    case {'datetime'}
        T = t;
        % T.TimeZone = TZ;
    case {'struct'}
        % t.TimeZone = TZ;
        UTCoffset = hours(tzoffset(t));
        D = datevec(t);
        T = struct('year',D(:,1),'month',D(:,2),'day',D(:,3),...
                   'hour',D(:,4),'minute',D(:,5),'second',D(:,6),...
                   'UTCOffset',UTCoffset);
    case {'num','datenum'}, T = datenum(t);
    case {'vec','datevec'}, T = datevec(t);
    otherwise
        error('Unrecognized output format');
    end

    if ~isequal(size(T),sz) && numel(T) == prod(sz), T = reshape(T,sz); end
end

function T = any2datetime(X,fmt,varargin)
% Convert X to datetime D, change X.TimeZone to TZ or X.UTCOffset, if available.
    
    if isa(X,'datetime'), T = X; return; end 
    
    if isempty(X), T = datetime.empty; return; end
    
    if isstruct(X), T = struct2datetime(X); return; end

    if ischar(X), X = cellstr(X); end
    if iscellstr(X) || isstring(X)
        if isempty(fmt)
            fmt = {};
            % avoid UnsupportedSymbol error on ISO 8601 date-stamps
            if contains(X{1},'T'), X = strrep(X,'T',' '); end 
        else
            fmt = {fmt};
        end
        T = datetime(datenum(X,fmt{:}),'convertfrom','datenum');
        return;
    end
    
    assert(isnumeric(X),'Unrecognized time format')
    
    if ~isempty(fmt)
        switch fmt
        case {'datetime','struct'}
            error('Input does not match requested convertion from %s',fmt); 
        case {'vec','datevec'}
            assert(ismatrix(X) && any(size(X,2) == [3,6]) && ...
                ~any(mod(X(:,1:min(size(X,2),5)),1) > 0,'all'),'Input is not valid DATEVEC');
            T = datetime(X);
            return;
        case 'num', fmt = 'datenum';
        end
    end
    if ~isempty(fmt)
        T = datetime(X,'convertfrom',fmt,varargin{:});
        return;
    end
    
    most = @(x) logical(mode(x));
    
    if ismatrix(X) && any(size(X,2) == [3,6]) && ~any(mod(X(:,1:min(size(X,2),5)),1) > 0,'all')
        T = datetime(X);
    elseif most(X > 1e9)
        warning('parsetime:mposix','Assuming numeric values are 1000·POSIX time');
        T = datetime(X/1000,'convertfrom','posixtime');
    elseif most(X > 1e6)
        warning('parsetime:posix','Assuming numeric values are POSIX time');
        T = datetime(X,'convertfrom','posixtime');
    elseif most(X < 1e5)
        warning('parsetime:excel','Assuming numeric values are EXCEL time');
        T = datetime(X,'convertfrom','excel');
    else
        T = datetime(X,'convertfrom','datenum');
    end
end

function t = struct2datetime(X)
% Convert time-structure (e.g. output of PVL_MAKETIMESTRUCT) into DATEVEC @ UTC, storing
% TimeZone or UTCOffset field

    assert(isscalar(X),'PARSETIME doesn''t work with non-scalar structures');
    if isfield(X,'UTCOffset')
        assert(isnumeric(X.UTCOffset),'Expecting numeric X.UTCOffset. Did you mean X.TimeZone?');
        UTCoffset = X.UTCOffset;
        X = rmfield(X,'UTCOffset');
        if ~isfield(X,'TimeZone')
            u = mode(UTCoffset);
            X.TimeZone = sprintf('%+03d:%02d',fix(u),floor(mod(u,1)*60)); % '+##:##'
            UTCoffset = UTCoffset - u;
        elseif ~any(strcmpi(X.TimeZone,{'UTC','+00:00','-00:00'}))
            warning('Applying both UTCOffset AND TimeZone correction');
        end
    else
        UTCoffset = 0;
    end
    if isfield(X,'DoY'), X = rmfield(X,'DoY'); end
        
    DEF = cell2struct([num2cell(datevec(floor(now))) {'UTC'} ]',...
        {'year','month','day','hour','minute','second','TimeZone'});
    X = completestruct(X,DEF,'warning','B~A');

    t = datetime(X.year,X.month,X.day,X.hour,X.minute,X.second,'TimeZone',X.TimeZone);

    if any(UTCoffset)
        t = t - hours(UTCoffset);
    end
end

function test()

    t = datetime(1985,03,10,23,30,0) + minutes(0:5:60)'; % most important hour in history!
    t.TimeZone = 'America/Mexico_City';
    
    % Basic run
    [T,dt,idx] = parsetime(t);
    assert(isequal(t,T) && isequal(dt,minutes(5)) && isequal(idx,(1:13)'));
    
    % Introduce some noise, remove points
    r = t + minutes(1.5 + rand(size(t))); 
    r(5:6) = [];
    [T,dt,idx] = parsetime(r,'-regular','latticeunit',minutes(5));
    assert(isequal(t,T) && isequal(dt,minutes(5)) && isequal(idx,([1:4,7:13])'));
    
    warning_disabler = naptime('cmd:reptimestamps'); %#ok<NASGU>
    
    % Test -unique and -sorted
    s = randi(numel(t),numel(t),1);
    [T,~,ic,ia] = parsetime(t(s),'-unique','-sorted');
    assert(isequal(T(ic(ia)),t(s)));
    
    % Convert back and forth
    types = {'datevec','datenum','struct'};
    for j = 1:numel(types)
        c = parsetime(t,'output',types{j});
        if isnumeric(c)
            r = parsetime(c,'convertfrom',types{j},'TimeZone',t.TimeZone);
            %r = parsetime(c);
        else
            r = parsetime(c); 
        end
        assert(isequal(t,r),'Test failed for %s',types{j});
    end
    
    fprintf('PARSETIME tests passed with flying colors\n');
end