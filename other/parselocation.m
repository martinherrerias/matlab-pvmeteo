function Loc = parselocation(Loc,varargin)
% LOC = PARSELOCATION(LOC) - make sure LOC is a structure with valid fields:
%
%   LOC.latitude (scalar, degrees north)
%   LOC.longitude (scalar, degrees east)
%   LOC.altitude (scalar, meters above sea level, optional)
%   LOC.name (string, optional)
%   LOC.TimeZone (string, passed to CHECKTIMEZONE, optional)
%
% LOC = PARSELOCATION(LOC,'optional',OPT) - control which fields are optional (just omitted if
%   missing/invalid). Default is {'TimeZone','altitude','name'}. Any missing/invalid field that
%   is not optional will cause an error. 'optional','all' makes all fields optional.
%
% LOC = PARSELOCATION(LOC,'-useweb') - use GEONAMES(LOC.latitude,LOC.longitude) to retrieve
%   altitude, TimeZone, and/or name (nearest toponym) and complete missing/invalid fields.
%
% LOC = PARSELOCATION(LOC,'-soft') - uses TimeZone = round(LOC.longitude/15), altitude = 0,
%   and/or name = 'unknown' (with warnings) instead of throwing an error. It is different from
%   using 'optional','all' in that the resulting structure is guaranteed to have all fields.
%
% '-inland' - issue a warning if the provided origin is not inside land polygons (10 km offset and
%   simplified from Matlab's 'coastlines').
%
% '-plot' - plot point and land polygons in Mercator projection.
%
% TODO: turn into a class, parse upon construction and pass around without having to parse again
%
% See also: PARSETIME, GEONAMES

    ALIAS = {'latitude',{'lat','latitude'};
            'longitude',{'lon','longitude'};
            'altitude',{'alt','altitude'};
            'TimeZone',{'timezone','time-zone','utcoffset','utc-offset'};
            'name',{'name','id'}};
    KNOWN = ALIAS(:,1);
    
    Loc = renamefields(Loc,ALIAS(:,[2,1]),'-ignorecase');
    if ~all(isfield(Loc,{'latitude','longitude'}))
    % Try using other fields, like origin, system, etc.
        try
            [Loc.longitude,Loc.latitude] = checkcoordsystem(0,0,0,Loc,'output','abs');
        end
    end
    
    [opt,varargin] = getflagoptions(varargin,{'-useweb','-soft','-inland','-plot'});
    opt.optional = {'TimeZone','altitude','name'};
    opt = getpairedoptions(varargin,opt,'restchk');
    opt.optional = parselist(opt.optional,ALIAS,'option');

    if opt.soft
        required = {};
    else
        required = setdiff(KNOWN,opt.optional);
    end
    
    if isempty(Loc) && ~isstruct(Loc), Loc = struct(); end

    if opt.useweb
    % Parse with all-optional fields, removing invalid
        Loc = parselocation(Loc,'optional',KNOWN,'useweb',false);
        
        missing = setdiff(KNOWN,fieldnames(Loc));
        if ~isempty(missing)
            try
                L = geonames(Loc.latitude,Loc.longitude);
                Loc = completestruct(Loc,L);
                warning('Using location %s from GEONAMES',shortliststr(missing))
            catch ERR2
                warning(ERR2.identifier,'Failed to use GEONAMES: %s',ERR2.message);
            end
        end
    end
   
    isnice = @(x) isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x);
    testfield('latitude',@(x) isnice(x) && abs(x) <= 90);
    testfield('longitude',@(x) isnice(x) && abs(x) <= 180);
    testfield('altitude',@(x) isnice(x) && x > -500 && x < 9000);
    testfield('TimeZone',@(x) ~isempty(checktimezone(x)));
    testfield('name',@ischar);

    if isfield(Loc,'TimeZone'), Loc.TimeZone = checktimezone(Loc.TimeZone); end
        
    if opt.soft
        if ~isfield(Loc,'TimeZone')
            Loc.TimeZone = checktimezone(round(Loc.longitude/15));
            warning('parselocation:timezone','Assuming location time zone UTC%s',Loc.TimeZone);
        end
        if ~isfield(Loc,'altitude')
            Loc.altitude = 0;
            warning('parselocation:altitude','Assuming altitude = 0 (sea level pressure)');
        end
    end
    
    if opt.inland || opt.plot
        L = load('landpolygons.mat','land');

        if opt.plot
            GUIfigure('map','Project Site',[0.1 0.1 0.6 0.6]); 
            clf(); axis equal; axis([-180 195 -80 90]);  hold on; 
            xticks(-165:15:180); yticks(-75:15:75); grid on;
            polyplot(polygon([-180 195],[-90 90]),[0.4 0.6 1 0.5],'none');
            polyplot(L.land,[0.4 0.8 0 0.5],[0 0.6 0.2]);
            plot(Loc.longitude,Loc.latitude,'rx');
        end

        if opt.inland && ~insidepolygon(L.land,Loc.longitude,Loc.latitude),...
            warningOrInfo('inland',...
                'Pretty sure your site is not on land. Verify coordinate system settings');
        end
    end

    function testfield(f,fun)
        isrequired = ismember(f,required);
        isthere = isfield(Loc,f);
        if ~isthere
            assert(~isrequired,'Missing location field: %s',f);
        else
            try 
                valid = fun(Loc.(f));
                assert(valid,'Invalid location field: %s',f);
            catch ERR
                valid = false; 
            end
            if ~valid
                if isrequired, throwAsCaller(ERR); end
                warning('Removing invalid location field: %s',f);
                Loc = rmfield(Loc,f);
            end
        end
    end
end