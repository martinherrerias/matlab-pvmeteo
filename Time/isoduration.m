function S = isoduration(dt,usemonths)
% S = ISODURATION(DT) - represent DURATION DT in ISO 8601 format, round to the nearest second.
% S = ISODURATION(DT,USEMONTHS) - By default, the function will only return months as a unit
%   (e.g. P1Y6M) if DT is an exact multiple of years(1)/12. That is, a month* is defined within
%   this function as 1/12th of a year (30.4369 days). Setting USEMONTHS = 1 will use this unit 
%   for any DT > 1 month*, while setting USEMONTHS = 0 will never use it.
%    
% DT = ISODURATION(S) - convert ISO 8601 string into DURATION scalar.

    narginchk(1,2); 
    if nargin < 2, usemonths = []; end
    
    if isa(dt,'duration') && isscalar(dt)
    	S = duration2string(dt,usemonths);
    elseif ischar(dt) || isstring(dt)
        S = string2duration(dt);
    else
        error('Unrecognized argument(s)');
    end
end

function S = duration2string(dt,usemonths)
   if isnan(dt), S = 'NA'; return; end
    if isempty(dt), S = duration.empty(); return; end

    dt = seconds(dt);
    if dt < 0.5, S = 'PT0S'; return; end
        
    S = 'P';
    x = floor(dt/31556952); % years
    if x > 0
        S = [S, num2str(x), 'Y']; 
        dt = dt - x*31556952;
    end
    if isempty(usemonths), usemonths = mod(dt,2629746) == 0; end
    if usemonths
        x = floor(dt/2629746); % months*
        if x > 0
            S = [S, num2str(x), 'M']; 
            dt = dt - x*2629746;
        end
    end
    x = floor(dt/86400); % days
    if x > 0
        S = [S, num2str(x), 'D']; 
        dt = dt - x*86400;
    end
    if dt == 0, return; end
    
    S = [S,'T'];
    x = floor(dt/3600); % hours
    if x > 0
        S = [S, num2str(x), 'H']; 
        dt = dt - x*3600;
    end
    x = floor(dt/60); % minutes
    if x > 0
        S = [S, num2str(x), 'M']; 
        dt = dt - x*60;
    end
    x = round(dt);
    if x > 0
        S = [S, num2str(x), 'S']; 
    end 
end

function DT = string2duration(S)
    S = upper(S);
    assert(numel(S) > 2,'Expecting character string with 3 or more elements');
    
    tokens = regexp(S,'^P(\d+Y)?(?<!.*T.*)(\d+M)?(\d+D)?T?(\d+H)?(\d+M)?(\d+S)?$','tokens');
    tokens = [tokens{:}];
    assert(~isempty(tokens),'Failed to match %s to ISO 8601 duration pattern',S);
    
    UNITS = hours([8765.82, 8765.82/12, 24, 1, 1/60, 1/3600]);
    used = ~cellfun(@isempty,tokens);
    
    DT = cellfun(@(x) str2double(x(1:end-1)),tokens(used))*UNITS(used)';     
end