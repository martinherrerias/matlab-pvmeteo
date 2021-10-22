function MD = filter(MD,varargin)
% MD = FILTER(MD,varargin)
%   Filtering & Downsampling of MeteoData object, as parsed by COMPLETEMETEODATA
%
% FILTER(..,'downsample',m) will reduce the time-step-resolution of the data
%       by a factor m, i.e. average data every m time-steps.
% FILTER(..,'downsample',DT) for duration DT will reduce the time-step-resolution of 
%       the data by m = DT/MD.timestep.
% 
%   NOTE: downsampling currently requires an integer m >= 1 (calls AVGDOWNSAMPLE)
%
% TODO: use (..,'resample',varargin) as wrapper for RESAMPLESTRUCTURE on kt/kd space
% TODO: use kt/kd space for matrix filters!
%
% FILTER(..,'filter',EXPR) will use boolean/numeric-index, keyword, or MATLAB-string-
%       expresion EXPR to select a sub-set of the data, and return filtered structures.
%       Numeric/boolean indices will work as long as X(EXPR) is valid for a data vector X.
%
%       Alternatively, EXPR can be a string-expression that works with eval(EXPR), using
%       FIELDNAMES(MD) or properties like t, dark, available, etc. e.g. 'GHI > 0' filters-out
%       all dark points; 'day(t) == 21 & month(t) == 3' reduces the data set to points for 
%       the 21st of march, and 'sunel > 10' selects points with solar elevation > 10°.
%
%       Finally, following are recognized keywords:
%           'all' - uses the complete time-series - equivalent to [], '1:Nt', or 'true(Nt,1)'
%           'notdark' - equivalent to 'dark ~= 0' 
%           'wholedays' - returns also the timesteps immediately before and after any 'notdark'
%           'Xdays' - select X uniformly spaced whole days
%           'Xrnd' - select X random samples, e.g. '100rnd' for a small test.
%
% See also: COMPLETEMETEODATA

    [opt,varargin] = getflagoptions(varargin,{'-verbose'});
    opt.filter = [];
    opt.downsample = [];
    opt.minavail = 0.5;
    [opt,varargin] = getpairedoptions(varargin,opt);

    if isempty(opt.filter) && isempty(opt.downsample)
       warning('Nothing to do!');
       return;
    end
    if isempty(opt.downsample), opt.downsample = 1; end

    printif = @(varargin) opt.verbose && fprintf(varargin{:});

    original_interval = MD.interval;
    if ~isregular(MD.data)
        MD = MeteoData.loadobj(MD);
    end
    if ~any(MD.interval=='ci'), original_interval = MD.interval; MD.interval = 'c'; end
    dt = MD.timestep;

    localtime = MD.t;
    localtime.TimeZone = MD.location.TimeZone;

    usefulsteps = parsefilter(MD,opt.filter);

    % Average data every 'timestep' points to get desired resolution
    if isduration(opt.downsample), opt.downsample = opt.downsample/dt; end
    if opt.downsample > 1

        printif('\nDownsampling x%d to %0.1f min.\n',opt.downsample,minutes(dt).*opt.downsample);

        available = MD.available;

        MD.t.TimeZone = 'UTC';
        MD.interval = 'c';
        
        MD.data = avgdownsample(MD.data,opt.downsample,usefulsteps,'logical',opt.minavail,'omitnan',varargin{:});
        if ~isempty(MD.uncertainty)
           MD.uncertainty = avgdownsample(MD.uncertainty.data,usefulsteps,'omitnan');
        end
        MD.flags = avgdownsample(MD.flags,opt.downsample,usefulsteps);
        usefulsteps = avgdownsample(usefulsteps,opt.downsample,'logical',opt.minavail,varargin{:});

        % Average Time structure (not considering filter)
        % t = avgdownsample(MD.t,opt.downsample); % UTC

        MD.timestep = dt*opt.downsample;
        MD.info{end+1} = sprintf('%dX-Downsampled',opt.downsample);

        % Update effective solar position, using Linke-Turbidity data, if available
        % TODO: let solar position be irradiance weighted instead, i.e. try to keep closure
        fld = {'sunel','sunaz','hourangle'};
        fld = fld(isfield(MD,fld));
        if ~isempty(fld)
            MD = rmfield(MD,fld);
            MD = getsunpos(MD);
        end

        if MD.interval == 'i', MD.interval = 'c'; end
    elseif opt.downsample < 1
        warning('Cannot yet resample');
    end

    if ~all(usefulsteps)
    % Clean Up: eliminate all useless time steps from MD and Time structures 

        if opt.downsample > 1
            printif('Keeping %d %dX downsampled steps out of %d available data points\n',...
                nnz(usefulsteps),opt.downsample,nnz(available));
        else
            printif('Keeping %d of %d available data points\n',...
                nnz(usefulsteps),nnz(MD.available));
        end

        MD = filterstructure(MD,usefulsteps,varargin{:});
        MD.missing = ~usefulsteps; 
        % MD.t = t(usefulsteps);
    end

    % Bring time-labels back to the original summarization
    if ~any(original_interval=='ci'), MD.interval = original_interval; end
end

function filter = parsefilter(MD,expr)
% Filter time steps for partial simulation (if required)

    MINAVAIL = 0.0;
    GOODAVAIL = 0.75;

    if isvector(expr) && (islogical(expr) && numel(expr) == MD.Nt) || ...
       isvector(expr) && (isnumeric(expr) && max(expr) <= MD.Nt)
   
        filter = expr(:); return;
        
    elseif ismatrix(expr) && size(expr,2) == MD.Nt, filter = expr;  return;     
    elseif isempty(expr) || isequal(expr,'all'), filter = true(MD.Nt,1);  return;
    elseif isempty(MD), filter = []; return;
    end
    
    validateattributes(expr,{'char'},{'nonempty'},'','filter expression');

    Nt = MD.Nt;
    t = MD.t;
    
    % [~,~,dayidx] = unique([Time.year,Time.month,Time.day],'row');
    [~,~,dayidx] = unique([year(MD.t),month(MD.t),day(MD.t)],'row');
    
    availability = accumarray(dayidx,MD.available,'sum')./ ...
                   accumarray(dayidx,~MD.dark,'sum');
               
    available = availability(dayidx) > MINAVAIL;

    notdark = ~MD.dark;
    wholedays = diff(notdark)~=0;
    wholedays = [wholedays;0] | [0;wholedays] | notdark;
    
    wholedays = wholedays & availability(dayidx) > GOODAVAIL;


    % Split keys of the form 'NNXX' into filter = 'XX', n = NN
    fkey = regexp(lower(expr),'(\d*)(\D+)','tokens');
    [n,fkey] = deal(fkey{1}{:});
    n = str2double(n);

    switch fkey
        case 'available', filter = available;
        case 'wholedays', filter = wholedays;
        case 'days'
        % XXdays - pick n uniformly spaced whole days
            assert(n > 0 && mod(n,1) == 0,'NNdays filter requires non-zero integer NN');
            f = MD.available;
            ic = zeros(Nt,1);
            [~,~,ic(f)] = unique([year(t(f)),month(t(f)),day(t(f))],'rows');
            if max(ic) < n
                warning('cmd:Ndays','No %d days are available, returning all %d',n,max(ic));
                filter = wholedays;
            else
                n = round(((1/n:1/n:1) - 0.5/n)*max(ic));
                filter = wholedays & any(ic == n,2);
            end
        case 'rnd'
        % XXrnd - pick n (predictable!) random samples
            assert(n > 0 && mod(n,1) == 0,'NNrnd filter requires non-zero integer NN');
            Na = nnz(MD.available);
            if Na < n
                warning('cmd:Nrnd','No %d samples are available, returning all %d',n,Na);
                n = Na;
            end
            f = false(Na,1);
            rng(12345); % make the sequence predictable, given Na,n
            f(sort(randperm(Na,n))) = true;
            filter = false(Nt,1);
            filter(MD.available) = f;
        otherwise
        % try to run whatever is in filter as matlab code
        
            pattern = cellfun(@(f) ['(?<!\.)(' f ')'],fieldnames(MD),'unif',0);
            expr = regexprep(expr,pattern,{'MD.data.$1'});
        
            try
                filter = eval(expr);
            catch ERR
                error('Tried and failed to evaluate filter = ''%s''. Got error: %s\n',...
                    expr,getReport(ERR));
            end
    end
end