function [loc,raw] = geonames(varargin)
% A simple interface for https://www.geonames.org.
% Requires a (free) user-account (see USERACCOUNTS.txt)
%
% [LOC,RAW] = GEONAMES('STR') - look for a place by name, and return the closest match. Makes a
%   web request of the form 'http://api.geonames.org/searchJSON?q=STR...' with additional custom
%   options (see below).
%
% LOC = GEONAMES(LAT,LON) - finds nearby toponym for coordinates LAT, LON (decimal degrees N, E).
%   Makes a web request of the form 'http://api.geonames.org/findNearbyJSON?lat=LAT&lng=LON&...'
%   with additional custom options (below).
%
% [..] = GEONAMES(..,'name',val) - provide additional search filters or custom options recognized 
%   by either web service, e.g. 'country','DE',... 'cities','cities15000',... etc. Values can be
%   passed as numbers, logicals or strings.
%
%   The following options are sets as defaults:
%
%   'maxRows',1 - By default GEONAMES returns only the closest match
%   'featureClass','P' - By default search is restricted to centers of population
%   'isNameRequired',true - By default search STR is set to be part of the places name.
%   'style','FULL' - return all information for the place including TimeZone and elevation data.
%
% OUTPUT:
%   LOC - struct (or array of struct) with fields {name,latitude,longitude,country} and optionally
%       (if style=FULL) {altitude, and TimeZone}.
%   RAW - struct (or array of struct) with all fields returned by the web service, including
%       info like administrative codes, population, etc. 
%
% REF: https://www.geonames.org/export/geonames-search.html
%      https://www.geonames.org/export/web-services.html
%
% EXAMPLES: 
%   [loc,raw] = geonames(48.7397,9.0971)                % HLRS
%   [loc,raw] = geonames('Leon, Guanajuato')            % my home town
%   [loc,raw] = geonames('Everest','featureClass','T')  % Mount Everest
%   [loc,raw] = geonames('Paris','maxRows',2)           % Paris, France; and Paris Texas

    COMMON = {'username','maxRows','featureClass','featureCode','cities','style'};
    SEARCH = {'q','name','name_equals','name_startsWith','startRow','country','countryBias',...
             'continentCode','adminCode1','adminCode2','adminCode3','adminCode4','adminCode5',...
             'tag','operator','charset','fuzzy','east','west','north','south',...
             'searchlang','orderby','inclBbox','isNameRequired','lang'};
    NEARBY = {'lat','lng','radius','localCountry'};
    ALL = unique([COMMON,SEARCH,NEARBY],'stable');

    DEF.maxRows = 1;
    DEF.username = '';
    DEF.featureClass = 'P';
    DEF.isNameRequired = true;
    DEF.style = 'FULL'; % SHORT,MEDIUM,LONG,FULL

    % Parse options without defaults
    [OPT,varargin] = getpairedoptions(varargin,ALL);
    OPT = completestruct(OPT,DEF);
    
    if isempty(OPT.username)
        try
            OPT.username = useraccounts('geonames');
            validateattributes(OPT.username,{'char'},{'nonempty'});
        catch
            error('Unable to retrieve valid user account');
        end
    end

    % this function was meant to filter out incompatible options, but the web server doesn't seem
    % to mind them, so just forward everything
    reduced = @(OPT,~) OPT;
    % reduced = @(OPT,fields) rmfield(OPT,setdiff(fieldnames(OPT),[fields,COMMON]));

    if numel(varargin) == 1 && ischar(varargin{1})
        url = 'http://api.geonames.org/searchJSON?';
        OPT.q = varargin{1};

        r = makerequest(url,reduced(OPT,SEARCH));
    elseif numel(varargin) == 2 && all(cellfun(@isnumeric,varargin))
        url = 'http://api.geonames.org/findNearbyJSON?';
        OPT.lat = varargin{1};
        OPT.lng = varargin{2};
        OPT = rmfield(OPT,'isNameRequired');

        r = makerequest(url,reduced(OPT,NEARBY));
    else
        error('Unrecognized syntax');
    end

    if ~isstruct(r) || ~isfield(r,'geonames') || isempty(r.geonames)
        if nargout > 1
            warning('Unexpected/empty results from web request');
            raw = r; loc = struct(); return;
        else
            error('Unexpected/empty results from web request');
        end
    end
    raw = r.geonames;
    if iscell(raw), raw = structarray(raw); end
    
    loc(numel(raw)) = struct();
    [loc.name] = deal(raw.toponymName);
    x = cellfun(@str2double,{raw.lat},'unif',0);
    [loc.latitude] = deal(x{:});
    x = cellfun(@str2double,{raw.lng},'unif',0);
    [loc.longitude] = deal(x{:});
    if strcmp(OPT.style,'FULL')
        [loc.altitude] = deal(raw.astergdem);
        x = arrayfun(@(x) x.timezone.timeZoneId,raw,'unif',0);
        [loc.TimeZone] = deal(x{:});
        [loc.country] = deal(raw.countryName);
    end
end

function r = makerequest(url,opt)
% Concatenate a web request address of the form URL?name=val&...&name=val
    args = struct2cell(opt);
    notchr = ~cellfun(@ischar,args);
    args(notchr) = cellfun(@(x) mat2str(x(:)),args(notchr),'unif',0);
    args(notchr) = strrep(args(notchr),'[','');
    args(notchr) = strrep(args(notchr),']','');
    args(notchr) = strrep(args(notchr),';',',');
    
    args = cellfun(@(a,b) [a '=' b],fieldnames(opt),args,'unif',0);
    call = [url,strjoin(args,'&')];
    r = webread(call);
end
