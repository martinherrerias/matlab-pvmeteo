function MD = completemeteodata(MD)
% MD = COMPLETEMETEODATA(MD) - a collection of data-processing and quality-control routines that
%   leave a MeteoData object ready to be used for PV-simulation.
%
%   - Check for errors in UTC-offset
%   - Calculate effective solar-position, and dependent quantities (clear-sky, indices, ..)
%   - Complete clear-sky and weather information from McClear / MERRA (if MD.options.webservices)
%   - Fill gaps in ancilliary variables (Ta, Patm,..) using MSSA.
%   - Perform data quality control, to detect inconsistencies, outliers and data-errors. 
%   - Plot rd-kt scatter, series heat maps, and print general information about the data.
%
% COMPLETEMETEODATA uses directly the options:
%
%   MD.options.verbose - if false, will skip the last step, not producing plots or messages.
%                        warnings will still be issued.
%
% 	MD.options.useweb - if true, COMPLETEMETEODATA attempts to retrieve CAMS McClear and MERRA2
%       data from CAMS WPS services (or previously downloaded files). CAMS data can be used for  
%       quality-check of irradiance measurements, and MERRA2 weather variables (RH, Patm, Ta, vw)
%       are used as defaults (if the given variable is missing) or to fill gaps in locally measured
%       variables.
%
%   MD.options.fitclearsky - if true (EXPERIMENTAL) try to fit a smooth series of Linke-Turbidty 
%       values to the data, to correct bias in clear-sky model. See FITCLEARSKY.
%   MD.options.clearskycheck - if true (EXPERIMENTAL) flag points with BNI > CSBNI·k1, and with
%       DHI < k2·CSDHI where factors k1, k2 are estimated from a inter-quartile-range of binned 
%       auto-detected clear-sky samples. (See METEOQC.CLEARSKY_TESTS).
%
% TODO:   Seasonal UTC-offset detection and correction
%         Shaded-sensor detection and correction
%         Integration of RBE / bestimate
%         Intra-hour resampling, for analysis of high-freq. modeling uncertainty 
%               (MPPT, transient temperature model, ramps, etc.)
%         User-Interactive tool for inspection, filtering, and correction
%
% See also: GUIMETEO, GETMETEODATA, PVLMOD_SPA, GTI2COMPONENTS, SENSORUNCERTAINTY
%   CHECKUTCOFFSET, METEOQC.TEST, METEODATA.REFRESH, GETSUNPOS

    opt = MD.options;

    printif = @(varargin) opt.verbose && fprintf(varargin{:});
    
    printif('Basic quality checks...\n');
    if ~isfield(MD,'ENI')
        MD = addsource(MD,'ENI',solarposition.extraradiation(MD.t),'sp.extraradiation');
    end
    MD = meteoQC.test(MD,'-basic');   % sunpos-independent quality checks

    MD = removemissing(MD,{'GHI'}); 
    
    printif('Checking UTC offset...\n');
    MD = checkUTCoffset(MD);  % verify Time Zone

    % from here on work in UTC, center-of-interval, uniform timesteps
    MD.t.TimeZone = 'UTC';
    if ~any(MD.interval == 'ic'), MD.interval = 'c'; end

    printif('Standard meteoQC.test...\n');
    MD = meteoQC.test(MD,{'BSRN','CIE'});
    MD = refresh(MD);
        
    MD = removemissing(MD); 

    printif('\n%s\n',flagsummary(MD.flags,'all',1));

    % Complete data structure with McClear + MERRA2
    if opt.useweb
        if ~all(isfield(MD,MeteoData.VARS_MERRA2))
            MD = getMERRA2(MD);
        end
        if ~all(isfield(MD,MeteoData.VARS_MCCLEAR))
            MD = getMcClear(MD);
        end
    end

    % TODO: should probably be its own function
    printif('Gap filling...\n');
    
    % gap-filling of ancilliary variables
    [b,MD.flags] = flagbit(MD.flags,{'interp'});

    fld = {'Patm','Ta','RH','vw','albedo','soiling'}; % check MeteoData.varnames('ambient')
    res = [8,1,1,1,24*7,24*7];

    for j = 1:numel(fld)
        if ~isfield(MD,fld{j}), continue; end

        bad = MD.flags.data.(fld{j}) ~= 0 & ~MD.missing;
        if ~any(bad), continue; end

        MD.flags.data.(fld{j}) = bitset(MD.flags.data.(fld{j}),b(1),bad);
        MD.data.(fld{j})(bad) = NaN;
        MD.data.(fld{j}) = MeteoData.fillgaps(MD.data.(fld{j}),MD.timestep,res(j),...
            'varname',fld{j},'exclude',MD.missing);
    end

    if ~isfield(MD,'Patm')
        MD = addsource(MD,'Patm',repmat(pvl_alt2pres(MD.location.altitude),MD.Nt,1),'pvl_alt2pres');
        % MD.flags.data.Patm = bitset(MD.flags.data.Patm,b(2));
        warning('cmd:patm','No pressure data provided: estimating from std. atmosphere');
    end

    if ~isfield(MD,'tpw') && all(isfield(MD,{'RH','Ta'}))
        
        warningdisabler = naptime('cleanup:average'); %#ok<NASGU>
        MDc = cleanup(MD); % use a sing Ta & RH channels
        tpw = pvl_calcPwat(MDc.data.Ta,MDc.data.RH*100)*10;
        clear MDc warningdisabler
        
        % if ~isfield(MD,'air_density'), rho = wetair.density(MD.Ta,MD.RH,MD.Patm); end
        MD = addsource(MD,'tpw',tpw,'pvl_calcPwat');
        % MD.flags.data.tpw = bitset(MD.flags.data.tpw,b(2));
    end

    % retrieve albedo & soiling from SimOptions, if missing
    fld = {'albedo','soiling'};
    for j = 1:numel(fld)
        if isfield(MD,fld{j}), continue; end
        [x,msg] = parsefactor(fld{j},MD.t); 
        MD = addsource(MD,(fld{j}),x,['SimOptions.',fld{j}]);
        warning([msg fld{j}]);
    end
    
    MD = getsunpos(MD);
    MD = refresh(MD);

    % Uncertainty-based bounds for closure
    MD = getuncertainty(MD);
    MD = bestimate(MD);
    MD = meteoQC.uncertainty_tests(MD);

    % EXPERIMENTAL: correct linke-turbidity from auto-detected clear-sky periods
    switch opt.fitclearsky
        case {1,'on'}, opt.fitclearsky = true;
        case {0,'off'}, opt.fitclearsky = false;
        case 'auto'
            opt.fitclearsky = all(contains(MD.getsourceof('CSGHI'),{'ineichen','esra','haurwitz'}));
        otherwise, error('Unknown option for fitclearsky');
    end
    if opt.fitclearsky
        MD = fitclearsky(MD,'verbose',opt.verbose,'minCSfraction',opt.minCSfraction);
    end
    if isfield(MD,'clearsky')
        MD.data.clearsky = MD.data.clearsky > 0;
        if opt.clearskycheck && nnz(MD.data.clearsky)/nnz(~MD.dark) > opt.minCSfraction
            MD = meteoQC.clearsky_tests(MD);
        end
    end

    % Plotting
    % if opt.verbose
    %     MD.t.TimeZone = MD.location.TimeZone;
    %     plot(MD,'ktrd','flagged',true);
    %     plot(MD,'heatmap');
    % end
end

function [MD,M2] = getMERRA2(MD)
% Complete data structure with MERRA2
% See: GETREMOTE

    M2 = struct();        
    % if ~MD.options.useweb
    %     warning('getMERRA2:disabled',...
    %         'Cannot retrieve MERRA2 data while obj.options.useweb is disabled'); 
    %     return; 
    % end
    if nargout < 2 && all(isfield(MD,MeteoData.VARS_MERRA2))
        warning('getMERRA2:allthere','All MERRA2 variables seem already available'); 
        return; 
    end
    printif = @(varargin) MD.options.verbose && fprintf(varargin{:});

    % Get MERRA2 1-hour weather data, resample to data resolution
    printif('Retrieving 1-hour MERRA2 data...\n');
%         try
        M2 = getremote('MERRA2',MD.location,MD.t,hours(1),'-sync',...
            'verbose',MD.options.verbose);
        M2.data = M2.data(:,ismember(M2.data.Properties.VariableNames,MeteoData.VARS_MERRA2));

        % u = M2.Properties.RowTimes;
        % assert(numel(MD.t) == numel(u) && all(abs(MD.t-u) < MD.timestep/2)); % DEBUG

        % MD.merra2.used = setdiff(M2.Properties.VariableNames,fieldnames(MD));
        % [M2,u] = timetable2struct(M2);

        M2.location = [];
        MD = [MD,M2];
        % MD = completestruct(MD,M2);

        printif('Successfully added to Meteo-Data\n');
%         catch ERR
%             warning(ERR.identifier,'Failed to retrieve MERRA-2 data! - %s',ERR.message);
%             printif('\n');
%         end
end
    
function [MD,MC] = getMcClear(MD)
% Complete data structure with McClear data
% See: GETREMOTE

    MC = struct();        
    % if ~MD.options.useweb
    %     warning('getMcClear:disabled',...
    %         'Cannot retrieve McClear data while obj.options.useweb is disabled'); 
    %     return; 
    % end
    if nargout < 2 && all(isfield(MD,MeteoData.VARS_MCCLEAR))
        warning('getMcClear:allthere','All McClear variables seem already available'); 
        return; 
    end
    printif = @(varargin) MD.options.verbose && fprintf(varargin{:});

    if MD.options.useweb && ~all(isfield(MD,MeteoData.VARS_MCCLEAR))
        % Get CAMS McClear 1-min Clear-Sky data for relevant period
        printif('Retrieving 1-min CAMS McClear data...\n');
        try
            MC = getremote('McClear',MD.location,MD.t,minutes(1),'-sync','verbose',MD.options.verbose);
            MC.data = MC.data(:,ismember(MC.data.Properties.VariableNames,MeteoData.VARS_MCCLEAR));

            MD = [MD,MC];

            % % MD.mcclear.used = setdiff(MC.Properties.VariableNames,fieldnames(MD));
            % [MC,u] = timetable2struct(MC);
            % assert(numel(MD.t) == numel(u) && all(abs(MD.t-u) < MD.timestep/2)); % DEBUG
            % MD = completestruct(MD,MC);

            printif('Successfully added to Meteo-Data\n');
        catch ERR
            warning(ERR.identifier,'Failed to retrieve McClear data! - %s',ERR.message);
            printif('\n');
        end
    end
end
 
function [v,msg] = parsefactor(x,t)
% V = PARSEFACTOR(X) - parse function-handle, scalar, or 12-vector X into a set of 
%   values for every time-step t. 
% V = PARSEFACTOR(KEY) - Do the same for X = getSimOption(KEY).
%
% Meant for seasonal/time-dependent factors like albedo, soiling, degradation, etc. 

    if ischar(x)
        [x,isdef] = getSimOption(x);
        if isdef, msg = 'Default '; else, msg = 'Custom '; end
    else
        msg = 'Custom ';
    end

    if isa(x,'function_handle')
        v = x(t);
        msg = [msg 'function ' func2str(x)];
    elseif isscalar(x)
        v = repmat(x,size(t)); 
        msg = sprintf('%sconstant (%0.2f)',msg,x); 
    elseif isnumeric(x) && numel(x)==12
    % Interpolate in time from the 12 monthly averages, using average days for Months,
    % as recommended from Klein, 1977 - Duffie Beckman 3rd Ed.
    % Repeat first and last two months, for smooth interpolation near edges.
            t0 = solarposition.mdoy([11:12,1:12,1:2]) + 365.*[-1,-1,zeros(1,12),1,1];
            x = x([11,12,1:12,1,2]);
           % doy = t(:) - datenum(year(t),0,0);
            doy = days(t - dateshift(t,'start','year'));
            v = interp1(t0,x,doy,'pchip');
            msg = [msg 'monthly values'];
    elseif isnumeric(x) && numel(x) == numel(t)
        v = x(:);
        msg = [msg 'vector of values'];
    else
        error('Scalar, 12-vector, Nt-vector, or function handle expected for %s',x); 
    end
end

