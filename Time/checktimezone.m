function [tz,key] = checktimezone(tz,default)
% TZ = CHECKTIMEZONE(X) - Check if X is a valid time-zone or UTCoffset. Valid inputs should be
%   IANA time  zones (e.g. 'Europe/Berlin','UTC',..); strings like 'UTC+3','+01:30',..;
%   scalar hour numbers (-11,3,10.5,..) and the keywords {'keep','none'} used by PARSETIME
%
% TZ = CHECKTIMEZONE(X,DEF) - Issue a warning if X is empty, and replace by (non-empty) DEF
%
% [TZ,KEY] = CHECKTIMEZONE(X,..) - return also char KEY (used by PARSETIME):
%
% X = 'keep' results in TZ = '', KEY = 'K' (keep)
% X = 'none' results in TZ = '', KEY = 'S' (set)
% X = '' results in TZ = DEF,    KEY = 'D' (default)
% otherwise parse TZ as usual,   KEY = 'S' (set)
%
% See also: PARSETIME, DATETIME, TIMEZONES

    if nargin < 2, default = ''; end
    if ~isempty(default), default = checktimezone(default); end

    key = 'S';
    if isempty(tz), tz = ''; key = 'D'; end
    
    if ischar(tz)
        tz = strtrim(tz);
        if isempty(tz)                 
            tz = default;
            key = 'D'; 
            if ~isempty(default)
                warning('checktimezone:default','Assuming time-zone: %s',default);
            end
            return; 
        end
        
        switch lower(tz)
            case 'keep', tz = ''; key = 'K'; return;
            case 'none', tz = ''; return;
            case {'z','utc','gmt','universal time','universal time (ut)'}, tz = 'UTC'; return;
        end

        % if contains(tz,'universal time','ignorecase',true)
        %     tz = 'UTC';
        %     return;
        % end

        % UTC+#:## / +##:## / +#.## ...
        s = regexp(tz,'^(UTC)?\s?(?<hour>[+\-\.0-9]{1,6}):?(?<min>\d{0,2})$','names');
        if ~isempty(s)
            if isempty(s.min), s.min = '0'; end
            tz = numerictimezone(str2double(s.hour) + str2double(s.min)/60);
            return;
        end
        
    elseif isa(tz,'duration'), tz = numerictimezone(hours(tz)); return;
    elseif isnumeric(tz), tz = numerictimezone(tz); return;
    else
        error('Expecting character array or offset hours (scalar)');
    end
    
    % Everything else, including IANA time zones
    try 
        t = datetime('now','TimeZone',tz);
        tz = t.TimeZone;
    catch ERR
        error('checktimezone:unknown','Unknown time zone');
    end
end

function tz = numerictimezone(k)
    assert(isscalar(k) && isreal(k) && isfinite(k) && ...
        mod(k,0.25) == 0 && k >= -12 && k <= 14,'checktimezone:range','Invalid time-zone');
    tz = sprintf('%+03d:%02d',fix(k),floor(mod(k,1)*60));
end