classdef (InferiorClasses = {?MeteoData}) meteoQC < table_ish
% A meteoQC object has an uint32 table column of flag values for each variable 
% (table column or structure field) of the object that is used for its construction. 
%
% See meteoQC_example.mlx 

properties (GetAccess=public,SetAccess=public)
    flags = {};
end
properties (Hidden,Constant)

    % TODO: use flag pairs instead? e.g. test+high / test+low
    % e.g. from https://wiki.pangaea.de/wiki/BSRN_Toolbox#Quality_Check
    %
    %  	1 	Measurement falls below the physically possible limit
    %  	2 	Measurement exceeds the physically possible limit
    %  	4 	Measurement falls below the extremely rare limit
    %  	8 	Measurement exceeds the extremely rare limit
    %  	16 	Compared with a corresponding measurement the measurement is too low
    %  	32 	Compared with a corresponding measurement the measurement is too high 
        
    FLAGS = {'NA','num','abs_phys','rel_phys','flat','interp','manual','other'...
             'BSRN_abs_lo','BSRN_abs_hi','BSRN_rare_lo','BSRN_rare_hi',...
             'CIE_lo','CIE_hi','CIE_out',...
             'CS_lo','CS_hi',...
             'UNC_hi','UNC_lo' };
    BITS = 32; % adaptive in the future?
    
    THRESH = 0.75;
end
methods (Static = true)
    function typ = TYPE(), typ = num2str(meteoQC.BITS,'uint%d'); end

    function MD = test(MD,varargin)
    % MD = METEOQC.TEST(MD,[SEQ,...])

        % const properties?
        ALLTESTS = {'NA','num','sunpos','weather','BSRN','CIE','clearsky','UNC'};
        FLD = ['ENI';MeteoData.varnames({'irradiance','ambient'})];
        BASIC = {'NA','num','sunpos','weather','BSRN','CIE'}; 
        MIN_FLAGS = {'NA','num','abs_phys'};
        
        OPT = {'sequence','basic','independent','iterations','valid','flag2nan'};

        validateattributes(MD,{'MeteoData'},{'nonempty'});

        % opt.independent = false;
        % opt.sequence = 'all';
        % opt.valid = @(x) isfinite(x);
        % opt.flag2nan = {'num'};
        % opt.iterations = 1;
        
        opt = getSimOption(cellfun(@(f) ['meteoQC.' f],OPT,'unif',0));
        opt = parseoptions(varargin,{'-independent','-basic'},opt.meteoQC,'dealrest');
        
        if isempty(MD.flags) || any(any(isfield(MD,FLD) & ~isfield(MD.flags,FLD)))
           MD.flags = meteoQC(MD,MD.flags); 
        end
        [~,MD.flags] = flagbit(MD.flags,MIN_FLAGS);
        
        opt.sequence = parselist(opt.sequence,ALLTESTS','test','in sequence');
        opt.flag2nan = parselist(opt.flag2nan,meteoQC.FLAGS','-soft');
        [opt.flag2nan,MD.flags] = flagbit(MD.flags,opt.flag2nan);
        
        if opt.independent, opt.iterations = 1; end
        
        if opt.basic
            if ~isempty(setdiff(opt.sequence,BASIC)) || ...
                    opt.iterations > 1
                % warning('-basic flag will disable other custom options')
                opt.iterations = 1;
                opt.sequence = intersect(opt.sequence,BASIC,'stable');
            end
            sunposargs = {{'ENI','TL'},0};
        else
            sunposargs = {};
        end

        for k = 1:opt.iterations
            for j = 1:numel(opt.sequence)
                switch opt.sequence{j}
                    case 'NA', MD.flags = checkfield(MD.flags,MD,@(x) ~isnan(x),'NA','all',true);
                    case 'num', MD.flags = checkfield(MD.flags,MD,opt.valid,'num','all',true);
                    case 'sunpos', MD = meteoQC.sunpos(MD,sunposargs{:});
                    case 'weather', MD = meteoQC.weather_tests(MD);
                    case 'BSRN', MD = meteoQC.BSRN_tests(MD);
                    case 'CIE', MD = meteoQC.CIE_tests(MD);
                    case 'UNC', MD = meteoQC.uncertainty_tests(MD);
                    case 'clearsky', MD = meteoQC.clearsky_tests(MD);
                    otherwise
                end
                
                if ~opt.independent && ~isempty(opt.flag2nan)
                    MD = meteoQC.flagged2nan(MD,opt.flag2nan); 
                end
            end
        end
        
        if opt.independent && ~isempty(opt.flag2nan)
            MD = meteoQC.flagged2nan(MD,opt.flag2nan); 
        end
    end
    
    function MD = sunpos(MD,fld,full)
    % MD = SUNPOS(MD,[FLD,FULL]) - Check solar-position, and solar-position-dependent variables
   
        FLD = {'ENI','TL','sunel','sunaz','hourangle','declination'};
        
        validateattributes(MD,{'MeteoData'},{'nonempty'});
        
        % if isfield(MD,'sza') && ~isfield(MD,'sunel')
        %    MD.sza = 90-MD.sza;
        %    MD = renamevars(MD,'sza','sunel');
        % end
        
        if nargin < 2 || (~iscell(fld) && isempty(fld)) || isequal(fld,'all')
            fld = intersect(FLD,fieldnames(MD));
        else
            % [ignored,ia] = setdiff(fld,fieldnames(MD));
            % fld(ia) = [];
            % if any(ia)
            %    warning(shortliststr(ignored,'Ignoring missing field')); 
            % end
            fld = intersect(fld,fieldnames(MD));
        end
        if isempty(fld), return; end
        
        ERRW = 7.5; % up to 7.5° degrees error due to interval summarization in hourly intervals
        ERRD = 2.0; % up to 2° minimum error for declination, refraction, etc.
        NPTS = 500;
        MPTS = 20;
        
        if nargin < 3 || isempty(full), full = (MD.Nt <= 2*NPTS); end

        MD.flags = checkfield(MD.flags,MD,@(x) ~isnan(x),'NA',fld,true);
        [B,F] = flagbit(MD.flags,{'num','abs_phys','rel_phys','model'});
        
        if contains(fld,'ENI')
        % Allow both the old solar constant of 1366 W/m² as well as the lowest values from 
        % SORCE/TIM, with some margin for modeling error. Replace flagged values with own 
        % estimates.
            kr2 = solarposition.extraradiation(MD.t,1);
            ok = MD.ENI < kr2*1369 & MD.ENI > kr2*1355;
            F.data.ENI = meteoQC.check(MD.ENI,@(x) ok,B(4),F.data.ENI);
            MD.ENI(~ok) = kr2(~ok)*1361.6;
        end

        if contains(fld,'TL')
        % Flag Linke-Turbidity values < 2 as rare, and TL < 1 as unphysical
            F.data.TL = meteoQC.check(MD.TL,@(x) x > 1,B(2),F.data.TL);
            F.data.TL = meteoQC.check(MD.TL,@(x) x > 2,B(3),F.data.TL);
        end
        
        MD.flags = F;
        
        fld = setdiff(fld,{'ENI','TL'});
        if isempty(fld), return; end

        idx = find(~MD.missing);
        if numel(idx) < NPTS, idx = (1:MD.Nt)'; end
 
        if full || numel(idx) > MPTS
        % Recalculate position [on a small sample], to check for angle units, conventions, etc.
        
            if ~full
                idx = randperm(numel(idx),min(NPTS,numel(idx)));
            end
            t = parsetime(MD.t(idx),'interval',{MD.interval,'c'},'step',MD.timestep);
            [R.sunaz, R.sunel,~,ST,R.declination] = pvlmod_ephemeris(t,MD.location,1e5,10);
            R.hourangle = (ST - 12)*15;
            
            % error in cos(zenith) ~ d(cos(z))/dw·ERRW
            dcz = abs(cosd(MD.location.latitude).*cosd(R.declination).*sind(R.hourangle));
            dcz = max(dcz,ERRD*pi/180);
            
            % error in cos(sunaz) ~ d(cos(az))/dz·(dz/dw)·ERRW
            dca = abs((sind(MD.location.latitude).*secd(R.declination) + ...
                       sind(R.sunel).*tand(R.declination)).*secd(R.sunel).^3).*dcz;
            dca = max(dca,ERRD*pi/180);
            
            anglediff = @(ab) abs(ab-360*round(ab/360));
            
            for j = 1:numel(fld)
                
                n = nnz(isfinite(MD.(fld{j})(idx)) & isfinite(R.(fld{j})));
                if n < MPTS, continue; end

                switch fld{j}
                   case 'hourangle'
                       units = {'degrees','negative-degrees','radians','negative-radians'};
                       convf = {1,-1,180/pi,-180/pi};
                       okfcn = @(x) anglediff(x - R.hourangle) < ERRW;
                   case 'declination'
                       units = {'degrees','negative-degrees','radians','negative-radians'};
                       convf = {1,-1,180/pi,-180/pi};
                       okfcn = @(x) anglediff(x - R.declination) < ERRD;
                   case 'sunel'
                       units = {'degrees','radians'};
                       convf = {1,180/pi};
                       okfcn = @(x) anglediff(sind(x) - sind(R.sunel)) < dcz;
                   case 'sunaz'
                       azconv = {'N2E','S2W','N2W','S2E','E2N','W2S'};
                       units = {};
                       convf = {};
                       for k = numel(azconv):-1:1
                           convf{k,1} = @(x) solarposition.fixazimuth(x,'N2E',azconv{k});
                           units{k,1} = ['deg ' azconv{k}];
                           convf{k,2} = @(x) solarposition.fixazimuth(x*180/pi,'N2E',azconv{k});
                           units{k,2} = ['rad ' azconv{k}];
                       end
                       okfcn = @(x) anglediff(cosd(x) - cosd(R.sunaz)) < dca;
                end
                [~,ok,u] = meteoQC.check_units(MD.(fld{j})(idx),okfcn,units,convf,[],n/numel(idx)*0.75);
                F.data.(fld{j})(idx(~ok)) = bitset(F.data.(fld{j})(idx(~ok)),B(2));
                
                if u == 0
                    warning('Field %s seems incompatible with PVL_EPHEMERIS',fld{j});
                elseif u > 0
                % Apply unit conversion to full field
                    fcn = convf{u};
                    if isnumeric(fcn)
                        MD.(fld{j}) = fcn.*MD.(fld{j});
                    else
                        MD.(fld{j}) = fcn(MD.(fld{j}));
                    end
                end
            end
        else
            % warning?
        end
            
        if isfield(MD,'hourangle')
            MD.hourangle = solarposition.fixtoplusminus180(MD.hourangle);
        end

        if isfield(MD,'declination')
            MD.flags.data.declination = ...
                meteoQC.check(MD.declination,@(x) abs(x) <= 23.5,B(2),MD.flags.data.declination);
        end

        if isfield(MD,'sunel')
            MD.flags.data.sunel = ...
                meteoQC.check(MD.sunel,@(x) abs(x) <= 90,B(1),MD.flags.data.sunel);
        end

        if isfield(MD,'sunaz')
            MD.sunaz = solarposition.fixtoplusminus180(MD.sunaz);
        end
    end
 
    function MD = BSRN_tests(MD,basic)
    % MD = BSRN_tests(MD) - BSRN recoomended QC tests (absolute irradiance upper and lower bounds)
    % MD = BSRN_tests(MD,basic) - if basic = 1, ignore solar elevation (e.g. before UTC check)
    %
    % Ref: https://bsrn.awi.de/fileadmin/user_upload/bsrn.awi.de/Publications/
    %      BSRN_recommended_QC_tests_V2.pdf
    
        if nargin < 2, basic = false; end
        validateattributes(MD,{'MeteoData'},{'nonempty'});
        validateattributes(basic,{'numeric','logical'},{'scalar','binary'});
        
        FLD = intersect(fieldnames(MD),{'GHI','DHI','BNI','GTI','USW'});
        if isempty(FLD), return; end

        if basic || ~isfield(MD,'sunel'), cz = 1; 
        else
            cz = max(0,sind(MD.sunel));
            cz(~isfinite(MD.sunel)) = 1;
        end

        [B,F] = flagbit(MD.flags,{'BSRN_abs_lo','BSRN_abs_hi','BSRN_rare_lo','BSRN_rare_hi'});

        F = checkfield(F,MD,@(x) x > -4.0,B(1),FLD);
        F = checkfield(F,MD,@(x) x > -2.0,B(3),FLD);

        for fld = {'GHI','CSGHI'}
            if ~isfield(MD,fld{1}), continue; else, f = fld{1}; end
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < 1.5*MD.ENI.*cz.^1.2 + 100,B(2),F.data.(f));
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < 1.2*MD.ENI.*cz.^1.2 + 50,B(4),F.data.(f));
        end

        for fld = {'DHI','CSDHI'}
            if ~isfield(MD,fld{1}), continue; else, f = fld{1}; end
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < 0.95*MD.ENI.*cz.^1.2 + 50,B(2),F.data.(f));
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < 0.75*MD.ENI.*cz.^1.2 + 30,B(4),F.data.(f));
        end

        for fld = {'BNI','CSBNI'} 
            if ~isfield(MD,fld{1}), continue; else, f = fld{1}; end
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < MD.ENI,B(2),F.data.(f));
            F.data.(f) = meteoQC.check(MD.(f),@(x) x < 0.95*MD.ENI.*cz.^0.2 + 10,B(4),F.data.(f));
        end

        if isfield(MD,'USW')
            F.data.USW = meteoQC.check(MD.USW,@(x) x < 1.2*MD.ENI.*cz.^1.2 + 50,B(2),F.data.USW);
            F.data.USW = meteoQC.check(MD.USW,@(x) x < MD.ENI.*cz.^1.2 + 50,B(4),F.data.USW);
        end

        % There is no explicit recommendation for GTI, so this is just GHI + BNI
        % TODO: consider incidence angle
        if isfield(MD,'GTI')
            F.data.GTI = meteoQC.check(MD.GTI,@(x) x < MD.ENI.*(0.95*cz.^1.2 + 1) + 50,B(2),F.data.GTI);
            F.data.GTI = meteoQC.check(MD.GTI,@(x) x < MD.ENI.*(0.75*cz + 0.95).*cz.^0.2 + 40,B(4),F.data.GTI);
        end
        
        MD.flags = F;
    end

    function MD = CIE_tests(MD,basic)
    % MD = CIE_TESTS(MD) - CIE recommended QC tests: closure GHI = FHI + cos(z)BNI within 15%,
    %                      DHI < GHI within 10%, flagged as CIE_hi / CIE_lo
    % MD = CIE_TESTS(MD,basic) - if basic = 1, ignore solar elevation (e.g. before UTC check)
    % 
    %  NOTE: Points with sunel < 4° or GHI* < 20 W/m² that fail the tests are flagged as CIE_out
    %
    %  NOTE(*): it's not clear what to do in the case of redundant sensors, e.g. which DHI + BNI to 
    %  compare with which GHI. Currently a sensor is flagged only if it fails the test for every
    %  combination of the other (~nan) sensors, and the threshold of 20 W/m² is based on the max.
    %  of any available GHI and DHI + BHI.
    %
    % Younes, S., R. Claywell, and T. Muneer. 2005. “Quality Control of Solar Radiation Data: 
    %  Present Status and Proposed New Approaches.” Energy 30 (9): 1533–49.

        if nargin < 2, basic = false; end
        validateattributes(MD,{'MeteoData'},{'nonempty'});
        validateattributes(basic,{'numeric','logical'},{'scalar','binary'});
        
        % FLD = intersect(fieldnames(MD),{'GHI','DHI','BNI','GTI','USW'});
        % if isempty(FLD), return; end
        
        if ~isfield(MD,'GHI'), return; end

        if basic || ~isfield(MD,'sunel'), cz = []; 
        else
            cz = max(0,sind(MD.sunel));
        end

        [B,F] = flagbit(MD.flags,{'CIE_lo','CIE_hi','CIE_out'});

        G0 = max(0,max(MD.GHI,[],2));
        if ~isempty(cz) && all(isfield(MD,{'BNI','DHI'}))
        % Get all combinations of BHI + DHI, arrange in dimensions {Nt,GHI,BNI,DHI}
            GHIc = permute(MD.BNI.*cz,[1,3,2]) + permute(MD.DHI,[1,3,4,2]);
            G0 = max(G0,max(GHIc,[],3:4));
        else
            GHIc = [];
        end

        if isempty(cz)
            inrange = G0 >= 20;
            % warning('Possibly applying CIE QC checks outside recommended conditions'); 
        else
        % Flag values only if sunel > 4 and GHI > 20 W/m²
            inrange = ~(MD.sunel < 4 | G0 < 20);
        end
        % F = checkfield(F,MD,@(x) inrange | isnan(x) | (x == 0),B(3),FLD);

        if isfield(MD,'DHI')
        % Check GHI - DHI > 10% GHI, flag GHI and DHI as lo/hi otherwise

            LO_global = MD.GHI < permute(MD.DHI,[1,3,2]) - G0*0.1;
            notlo = ~all(LO_global,3);
            F.data.GHI = meteoQC.check(MD.GHI,@(x) ~inrange | notlo,B(1),F.data.GHI);
            F.data.GHI = meteoQC.check(MD.GHI,@(x) inrange | notlo,B(3),F.data.GHI);

            nothi = permute(~all(LO_global,2),[1,3,2]);
            F.data.DHI = meteoQC.check(MD.DHI,@(x) ~inrange | nothi,B(2),F.data.DHI);
            F.data.DHI = meteoQC.check(MD.DHI,@(x) inrange | nothi,B(3),F.data.DHI);
        end

        if ~isempty(GHIc)
        % Check closure: | BNI cz + DHI - GHI | < 15% GHI, flag each as lo/hi otherwise

            na = isnan(MD.GHI) | isnan(GHIc);
            LO_global = MD.GHI < 0.85*GHIc | na; % {Nt,GHI,BNI,DHI}
            HI_global = MD.GHI > 1.15*GHIc | na;

            notlo = ~all(LO_global,3:4);
            nothi = ~all(HI_global,3:4);
            F.data.GHI = meteoQC.check(MD.GHI,@(x) ~inrange | notlo,B(1),F.data.GHI);
            F.data.GHI = meteoQC.check(MD.GHI,@(x) ~inrange | nothi,B(2),F.data.GHI);
            F.data.GHI = meteoQC.check(MD.GHI,@(x) inrange | (nothi & notlo),B(3),F.data.GHI);

            notlo = permute(~all(HI_global,2:3),[1,4,2,3]);
            nothi = permute(~all(LO_global,2:3),[1,4,2,3]);
            F.data.DHI = meteoQC.check(MD.DHI,@(x) ~inrange | notlo,B(1),F.data.DHI);
            F.data.DHI = meteoQC.check(MD.DHI,@(x) ~inrange | nothi,B(2),F.data.DHI);
            F.data.DHI = meteoQC.check(MD.DHI,@(x) inrange | (nothi & notlo),B(3),F.data.DHI);

            notlo = permute(~all(HI_global,[2,4]),[1,3,2,4]);
            nothi = permute(~all(LO_global,[2,4]),[1,3,2,4]);
            F.data.BNI = meteoQC.check(MD.BNI,@(x) ~inrange | notlo,B(1),F.data.BNI);
            F.data.BNI = meteoQC.check(MD.BNI,@(x) ~inrange | nothi,B(2),F.data.BNI);
            F.data.BNI = meteoQC.check(MD.BNI,@(x) inrange | (nothi & notlo),B(3),F.data.BNI);
        end
        MD.flags = F;
    end

    function MD = uncertainty_tests(MD)
    % EXPERIMENTAL: Uncertainty-based QC

        validateattributes(MD,{'MeteoData'},{'nonempty'});
        if isempty(MD.uncertainty), return; end
        
        [B,F] = flagbit(MD.flags,{'UNC_lo','UNC_hi'});

        [~,MD0] = bestimate(MD);

        fields = intersect(fieldnames(MD),fieldnames(MD.uncertainty));
        fields = intersect(fields,fieldnames(MD0));
        for j = 1:numel(fields)
            notlo = MD.(fields{j}) >= MD0.(fields{j}) - MD.uncertainty.(fields{j});
            F.data.(fields{j}) = meteoQC.check(MD.(fields{j}),@(x) notlo,B(1),F.data.(fields{j}));

            nothi = MD.(fields{j}) <= MD0.(fields{j}) + MD.uncertainty.(fields{j});
            F.data.(fields{j}) = meteoQC.check(MD.(fields{j}),@(x) nothi,B(2),F.data.(fields{j}));
        end
        MD.flags = F;
    end

    function MD = clearsky_tests(MD,plotting)
    % EXPERIMENTAL: Clear-Sky-model-based QC

        validateattributes(MD,{'MeteoData'},{'nonempty'});
        
        REQ = {'CSGHI','CSBNI','CSDHI','ENI','clearsky'};
        if nargin < 2, plotting = true; end
        [B,F] = flagbit(MD.flags,{'CS_lo','CS_hi'});

        missing = ~isfield(MD,REQ);
        if any(missing)
            warning('Skipping clearsky test. %s',shortliststr(REQ(missing),'Missing field'));
            return;
        end

        dark = MD.dark;

        if nnz(MD.clearsky(~dark))/nnz(~dark) < MD.options.minCSfraction
            warning('Not enough clear samples for clear-sky test');
            return;
        end

        % Use best available (inverse-variance-weighted) estimate for BNI up to now.
        [~,CG,CU] = bestimate(MD);

        if isfield(MD,'BNI')
        % Check that BNI < CSBNI
            bad = MD.BNI > MD.CSBNI;
            if any(bad,'all')
                % Get an estimate of the minimum CSkn - kn, based on "sunny" samples ...
                % "sunny" =  clear-sky | any point with BNI > CSBNI
                f = MD.clearsky | CG.BNI > MD.CSBNI & ~dark;  
                CSkn = MD.CSBNI./MD.ENI;
                kn = CG.BNI(f)./MD.ENI(f);
                sn = min(CU.BNI(f)./MD.ENI(f));
                [lb,~,out] = binnedIQR(CSkn(f),(CSkn(f)-kn),'-asym','minrange',sn,'xlim',[0 1],'-interp');
                kn_max = CSkn - lb(CSkn);
                bad = MD.BNI./MD.ENI > kn_max; % ... to use as min. threshold

                F.data.BNI = meteoQC.check(MD.BNI,@(x) ~bad,B(2),F.data.BNI);
                MD.clearsky(f) = MD.clearsky(f) & ~out;

                if plotting
                    f(f) = ~out;
                    GUIfigure('clearskycheck','Clear-Sky Check','2:1'); clf(); 
                    subplot(1,2,1); hold on;
                    CSkn = repmat(CSkn,1,size(MD.BNI,2));
                    scatter(CSkn(:),reshape(MD.BNI./MD.ENI,[],1),1,[1 1 1]*0.8);
                    densityplot(CSkn(f,1),kn(~out),2);
                    plot(lb.GridVectors{1},lb.GridVectors{1}-lb.Values,'r:')
                    kn = MD.BNI./MD.ENI;
                    scatter(CSkn(bad),kn(bad),2,'r');      
                    grid on; plot([0,1],[0,1],'k:'); 
                    axis equal; axis([0 0.8 0 0.8]);
                    xlabel('Clear-Sky k_{n,CS} = BNI_{CS}/ENI');
                    ylabel('Measured k_n = BNI_{CS}/ENI');
                end
            end
        end

        if all(isfield(MD,{'DHI','GHI'}))
            % Check that CSrd < rd
            CSrd = MD.CSDHI./MD.CSGHI;
            bad = any(MD.DHI < CSrd.*CG.GHI,2) | any(CG.DHI < MD.GHI.*CSrd,2);
            if any(bad)
            % Get the minimum fraction rd/CSrd based on clear-sky samples

                f = MD.clearsky;
                rd = CG.DHI(f)./CG.GHI(f);
                sd = hypot(rd.*CU.GHI(f),CU.DHI(f))./CG.GHI(f); % conservative: assumes no correlation
                % [lb,~,out] = binnedIQR(CSkn(f),(CSkn(f)-kn),'-asym','minrange',sn,'xlim',[0 1],'-interp');
                % [lb,~,out] = binnedIQR(CSrd(f),rd./sd./CSrd(f),'-asym','xlim',[0 1],'minrange',1,'-interp','-plot');
                [lb,~,out] = binnedIQR(CSrd(f),rd./CSrd(f),'-asym','xlim',[0 1],'minrange',min(sd),'-interp');
                f(f) = ~out;
                rd =[MD.DHI./CG.GHI,CG.DHI./MD.GHI];
                rd_min = lb(CSrd).*CSrd; % ... minimum expected CSrd
                bad = rd < rd_min; % ... to use as min. threshold

                F.data.DHI = meteoQC.check(MD.DHI,@(x) MD.DHI./CG.GHI > rd_min,B(1),F.data.DHI);
                F.data.GHI = meteoQC.check(MD.GHI,@(x) CG.DHI./MD.GHI > rd_min,B(2),F.data.GHI);

                if plotting
                    GUIfigure('clearskycheck'); 
                    subplot(1,2,2); hold on;
                    CSrd = repmat(CSrd,1,size(rd,2));
                    scatter(CSrd(:),rd(:),1,[1 1 1]*0.8);
                    densityplot(CSrd(f,1),CG.DHI(f)./CG.GHI(f),2);
                    plot(lb.GridVectors{1},lb.Values.*lb.GridVectors{1},'r:')
                    scatter(CSrd(bad),rd(bad),2,'r');
                    grid on; plot([0,1],[0,1],'k:');
                    axis equal; axis([0 0.8 0 0.8]);
                    xlabel('Clear-Sky r_{d,CS} = DHI_{CS}/GHI_{CS}');
                    ylabel('Measured r_d = DHI/GHI_o, DHI_o/GHI')
                end
            end
        end
        MD.flags = F;
    end
        
    function MD = weather_tests(MD)
        
        validateattributes(MD,{'MeteoData'},{'nonempty'});
        [b,F] = flagbit(MD.flags,{'abs_phys','num','IQR'});
        
        FLD = {'Patm','Ta','RH','tpw','vw','albedo','soiling'};
        if ~any(isfield(MD,FLD)), return; end

        OPTS = {'P',MD.options.outlier_P,'sequence',{'num','IQR'},'-independent',...
            'flaglist',{'foo','num','IQR'}};
        
        dt = MD.timestep;
        if isempty(dt), dt = 0; else, dt = hours(dt); end
        
        for j = 1:numel(FLD)
            if ~isfield(MD,FLD{j}), continue; end
            
            switch FLD{j}
            case 'Patm'
                if ~isempty(MD.location) && isfield(MD.location,'altitude')
                    z = MD.location.altitude;
                else
                    z = [];
                end
                [ok,MD.Patm] = meteoQC.pressure_check(MD.Patm,z);
                constant = (MD.Patm == MD.Patm(1,:));
                f = minQC(MD.Patm,OPTS{:},'window',24/dt,'valid',@(x) ok);

            case 'Ta'
                [ok,MD.Ta] = meteoQC.temperature_check(MD.Ta);
                f = minQC(MD.Ta,OPTS{:},'window',3/dt,'valid',@(x) ok);

            case 'RH'
                [ok,MD.RH] = meteoQC.humidity_check(MD.RH);
                f = minQC(MD.RH,OPTS{:},'window',3/dt,'valid',@(x) ok);

            case 'tpw'
                [ok,MD.tpw] = meteoQC.tpw_check(MD.tpw);
                f = minQC(MD.tpw,OPTS{:},'window',3/dt,'valid',@(x) ok);

            case 'vw'
                [ok,MD.vw] = meteoQC.windspeed_check(MD.vw);
                f = minQC(MD.vw,OPTS{:},'window',2/dt,'valid',@(x) ok);
                
            case 'albedo'
                [ok,MD.albedo] = meteoQC.albedo_check(MD.albedo);
                f = minQC(MD.albedo,OPTS{:},'window',4/dt,'valid',@(x) ok);
                
            case 'soiling'
                [MD.(FLD{j}),ok] = meteoQC.check_units(MD.(FLD{j}),[0,1],{'','%'},{1,1e-2});
                f = minQC(MD.(FLD{j}),OPTS{:},'window',4/dt,'valid',@(x) ok);
                    
            end
            
            for k = 1:3
                F.data.(FLD{j}) = bitset(F.data.(FLD{j}),b(k),bitget(f,k));
            end
        end
        MD.flags = F;
    end

    function [ok,p] = pressure_check(p,z)
    % ok = meteoQC.pressure_check(P,z) - Check that atmospheric pressure P [Pa] is within
    %   WMO scale recommendations and historical pressure records. Use elevation z (mASL) to
    %   account for site elevation.
    % [ok,P] = meteoQC.pressure_check(P,z) - check & correct for alternate units.
    
        ABS = [0.5E5,1.08E5];        % WM0 scale recommendations [Pa]
        REL = [870,1084.8]/1013.25;  % barometric pressure records (Wikipedia)
        
        narginchk(1,3);

        validateattributes(p,{'numeric'},{'real','2d'});

        if nargout > 1
            [p,ok] = meteoQC.check_units(p,ABS,{'Pa','hPa','kPa','psi','bar'},{1,1e2,1e3,6894.76,1e5});
        else
            ok = p >= ABS(1) & p <= ABS(2);
        end

        if nargin > 1 && ~isempty(z)
            validateattributes(z,{'numeric'},{'real','scalar','>',-430,'<',6400});
            p0 = pvl_alt2pres(z);
            ok = p/p0 >= REL(1) & p/p0 <= REL(2);
        end
    end      
    function [ok,Ta] = temperature_check(Ta)
    % ok = meteoQC.temperature_check(Ta) - Check that air temperature Ta [°C] is within
    %            historical records. 
    % [ok,Ta] = meteoQC.temperature_check(Ta) - check & correct for alternate units.
    
        ABS = [-90,60];

        validateattributes(Ta,{'numeric'},{'real','2d'});
        
        if nargout > 1
            Ta = meteoQC.check_units(Ta,[-20,45],{'°C','K','°F'},{ 1,@(x) x - 273.15,@(x) (x-32)/1.8 });
        end
        ok = Ta >= ABS(1) & Ta <= ABS(2);
    end
    function [ok,RH] = humidity_check(RH)
    % ok = meteoQC.humidity_check(RH) - Check that relative humidity is within [0,1.2]
    % [ok,RH] = meteoQC.humidity_check(RH) - check & correct for alternate units.
    
        ABS = [0,1.2]; % allow some supersaturation

        validateattributes(RH,{'numeric'},{'real','2d'}); 
        if nargout > 1
            [RH,ok] = meteoQC.check_units(RH,ABS,{'0-1','%'},{ 1,1e-2 });
        else
            ok = RH >= ABS(1) & RH <= ABS(2);
        end
    end
    function [ok,tpw] = tpw_check(tpw)
    % [ok,tpw] = meteoQC.tpw_check(tpw) - Check total precipitable water content [mm ~ kg/m²]
    
        ABS = [0,80]; % Records seem to reach ~70-80 kg/m², but no official max.

        validateattributes(tpw,{'numeric'},{'real','2d'}); 
        if nargout > 1
            [tpw,ok] = meteoQC.check_units(tpw,ABS,{'mm','cm'},{1,10},[2,1]);
        else
            ok = tpw >= ABS(1) & tpw <= ABS(2);
        end
    end
    function [ok,vw] = windspeed_check(vw)
    % [ok,vw] = meteoQC.windspeed_check(vw) - Check wind speed [m/s at module level]
    %   The most common "unit" confusion here would be to use wind speed measured at another
    %   height, e.g. 10 m. A wind shear factor of 0.46 ~ log(2/0.5)/log(10/0.5) corresponding
    %   to a roughness length of 0.5 m is used as an approximation.
    
        ABS = [0,12.5]; % Storm force on Beaufort scale ("considerable structural damage").
         
        validateattributes(vw,{'numeric'},{'real','2d'}); 
        if nargout > 1
            [vw,ok] = meteoQC.check_units(vw,[0,5],{'m/s','m/s @ 10 m'},{1,0.46},[1,2]); %(§)
        else
            ok = vw >= ABS(1) & vw <= ABS(2);
        end
    end 
    function [ok,rho] = albedo_check(rho)
    % ok = meteoQC.albedo_check(rho) - Check that albedo values are within [0.04,0.9]
    % [ok,rho] = meteoQC.albedo_check(rho) - check & correct for alternate units (%)
    
        ABS = [0.04,0.9]; % fresh asphalt to fresh snow 
        
        validateattributes(rho,{'numeric'},{'real','2d'}); 
        if nargout > 1
            [rho,ok] = meteoQC.check_units(rho,ABS,{'','%'},{ 1,1e-2 });
        else
            ok = rho >= ABS(1) & rho <= ABS(2);
        end
    end

    function flag = check(value,okfcn,bit,flag,flagifnan)
    % F = meteoQC.check(V,fcn,B,F) - if ~fcn(V) & isfinite(V), set bit B on F
    % F = meteoQC.check(V,fcn,KEY,F) - ... set bit corresponding to KEY on F
    % F = meteoQC.check(V,fcn,B,F,1) - set B  if ~fcn(V), even if ~isfinite(V)
    %
    % See also: meteoQC.checkfield

        if nargin < 5, flagifnan = false; end
        if nargin < 4 || isempty(flag) || ( isscalar(flag) && flag == 0 ) 
            flag = zeros(size(value),meteoQC.TYPE); 
        end

        notok = ~okfcn(value);
        if ~flagifnan
            notok = notok & ~isnan(value);
        end
        flag(notok) = bitset(flag(notok),bit);
    end

    function  [x,pass,best_units] = check_units(x,okfcn,units,convfcn,order,threshold)
    % [x,ok,idx] = meteoQC.check_units(X,OKFCN,UNITS,CONVFCN,ORDER,THRESH) - check the units of  
    %   array X, by searching for IDX = argmax(j) OKFCN( CONVFCN{j}(x) ), where CONVFCN is a 
    %   cell-array  of unit-conversion functions, and OKFCN is a function that returns true for 
    %   values of X within a 'reasonably expected' range. 
    %
    % NOTES:
    %   + OKFCN can also be specified as a range [MIN,MAX]
    %   + CONVFCN{j} must take X [UNITS{j}] and return X[UNITS{1}]. It can also be a scalar factor.
    %   + The search is interrupted as soon as the fraction of 'reasonable' values reaches
    %     THRESHOLD,(default meteoQC.THRESH, use THRESHOLD = 1 to search accross all units).
    %   + A warning will be issued if the selected unit is not the first, i.e. UNITS{1}
    %   + The order in which units are tested can be customized by ORDER 
    
        narginchk(4,6);
        if nargin < 5 || isempty(order), order = 1:numel(units); end
        if nargin < 6 || isempty(threshold), threshold = meteoQC.THRESH; end

        validateattributes(x,{'numeric'},{'real','2d'});
        if isnumeric(okfcn)
            lim = okfcn;
            validateattributes(lim,{'numeric'},{'real','vector','numel',2,'increasing'},'','limits');
            okfcn = @(x) x >= lim(1) & x <= lim(2);
        else
            validateattributes(okfcn,{'function_handle'},{'scalar'},'','okfcn'); 
        end
        units = cellstr(units);
        m = numel(units);
        validateattributes(convfcn,{'cell'},{'numel',m,},'','factors'); 
        
        validateattributes(threshold,{'numeric'},{'scalar','positive','<=',1},'','threshold');
        validateattributes(order,{'numeric'},{'vector','integer','positive','<=',m,'numel',m},'','order'); 
        
        available = isfinite(x);
        pass = available;
        if ~any(available), return; end
        
        n = nnz(available);
        best_score = 0;
        best_units = 0;
        best_x = x;

        current_score = NaN;
        for j = order(:)'
            if isnumeric(convfcn{j})
                validateattributes(convfcn{j},{'numeric'},{'real','scalar','finite','nonzero'},'','factors{j}'); 
                u = x.*convfcn{j};
            else
                validateattributes(convfcn{j},{'function_handle'},{'scalar'},'','factors{j}'); 
                u = convfcn{j}(x);
            end
            ok = okfcn(u);
            if j == 1
                current_score = nnz(ok)/n;
            end
            if nnz(ok)/n > best_score
                best_score = nnz(ok)/n;
                best_units = j;
                best_x = u;
                pass = ok & available;
                if best_score > threshold, break; end
            end
        end
        
        x(available) = best_x(available);
        if best_units > 1
            warning('meteoQC:units',['Wrong units?: detected [%s] (%0.0f%% pass) ',...
                'converted to [%s] (%0.0f%% pass)'],...
                units{best_units},current_score*100,units{1},best_score*100);
        end
    end
    
    function MD = flagged2nan(MD,flag_bits,fields)
    % D = flagged2nan(D,flag_bits,fld) - replace elements in each MD.(fld) for which any
    %     bitget(F.(fld),flag_bits) > 0  by NaN.

        validateattributes(MD,{'numeric','MeteoData'},{'nonempty'});
        if isempty(flag_bits), return; end
        
        if ischar(flag_bits) || iscell(flag_bits)
            try
                flag_bits = flagbit(MD.flags,flag_bits);
            catch
                [flag_bits,idx] = parselist(flag_bits,MD.flags.flags,'-soft');
                if any(idx == 0)
                    warning('%s not found',shortliststr(flag_bits(idx == 0),'Flag','quotes',''''));
                    flag_bits = flag_bits(idx > 0);
                end
                flag_bits = flagbit(MD.flags,flag_bits);
            end
        end
            
        validateattributes(flag_bits,{'numeric'},{'vector','integer','positive'});
        if size(flag_bits,1) > 1, flag_bits = flag_bits'; end
        m = numel(flag_bits);
        flag_bits = shiftdim(flag_bits,-1); % [1,1,m]
    
        if nargin < 3 || (~iscell(fields) && isempty(fields))
            fields = fieldnames(MD);
        else
            fields = intersect(fields,fieldnames(MD)); 
        end

        for j = 1:numel(fields)
            V = MD.(fields{j});
            B = repmat(flag_bits,size(V));
            F = repmat(MD.flags.data.(fields{j}),1,1,m);
            if ~isequal(size(B),size(F))
               warning('Inconsistent flags for field: %s',fields{j});
               continue;
            end
            notok = any(bitget(F,B),3);
            if any(notok,'all')
                MD.(fields{j})(notok) = NaN;
            end
        end
    end
end
methods    
    function F = meteoQC(D,F0)
    % F = meteoQC(D) - Create an object similar* to D, with uint32 flag values (set to 0)
    % F = meteoQC(D,F0) - Initialize flag values for D from those in meteoQC object F0

        if nargin == 0, F.data = table.empty; return; end

        if nargin < 2, F0 = []; end
        if ~isempty(F0)
            assert(isa(F0,'meteoQC'),'Expecting existing meteoQC object');
        end
        
        if isa(D,'meteoQC') 
            m = max(F.data{:,:},'all');
            if m > uint64(2)^numel(F.flags) - 1
                error('More non-zero bits than flags?');
            end
            F = D;
            for j = 1:size(F,2)
                F.data{:,j} = typecast(F.data{:,j},meteoQC.TYPE);
            end
        elseif isa(D,'MeteoData')
            F = meteoQC(D.data,D.flags);
            if nargin > 1, F = merge(F,F0); end
            return;
        else
            if isstruct(D) && isscalar(D)
                fld = fieldnames(D);
            elseif isa(D,'tabular')
                fld = D.Properties.VariableNames;
            else
                error('Expecting table or structure');
            end

            C = cell(1,numel(fld));
            for j = 1:numel(fld)
                C{j} = zeros(size(D.(fld{j})),meteoQC.TYPE);
            end
            F.data = table(C{:},'VariableNames',fld);
        end
        
        if ~nargin > 1, F = merge(F,F0); end
    end
    
%     function obj = subsasgn(obj,s,varargin), obj = builtin('subsasgn',obj,s,varargin{:}); end
    
%     function disp(obj)
%         if numel(obj) == 0, builtin('disp',obj); return; end
%         header = [ regexprep(mat2str(obj.size),{'[',']',' '},{'','','x'}),...
%             ' <a href = "matlab:help meteoQC">meteoQC</a>'];
%         if numel(obj) == 0, disp([header newline()]); return; end
%         
%         header = [header ' with ' nthings(numel(obj.flags),'flag','zero','no')];
%         % if ~isempty(obj.flags)
%         %     header = [header ': ' shortliststr(obj.flags,'','quotes','''')];
%         % end
%         disp([header newline()]);
%         disp(obj.data);
%     end
    
    function [b,F] = flagbit(F,flag_keys)
    % b = flagbit(F,keys) - return bit indices corresponding to keys. Throw an error if any
    %   of the keys is not already in the the list F.flags.
    % [b,F] = append any new keys to the end of the list F.flags.
    %
    % NOTE: the bit order is defined by F.flags, namely: F.flags{j} <-> j'th MSB
    
        flag_keys = cellstr(flag_keys);
        [ia,b] = ismember(flag_keys,F.flags);
        if all(ia), return; end
        
        if nargout < 2
            error([shortliststr(flag_keys(~ia),'Unknown flag key','quotes',''''),...
                  ' use syntax [b,F] = flagbit(F,key) to expand flag list']);
        else
            neu = numel(F.flags)+(1:nnz(~ia));
            if neu(end) > meteoQC.BITS
                error(['meteoQC currently supports up to %d flags. ',...
                       'After belittling the author for his lack of foresight, ',...
                       'consider merging some of your flags (mergeflags), ',...
                       'using an additional meteoQC object for the same data, ',...
                       'or changing the constant meteoQC.BITS'],meteoQC.BITS);
            end
            F.flags(1,neu) = flag_keys(~ia);
            b = [b(ia),neu];
        end
    end
    
    function n = flags2code(obj,c)
        c = cellstr(c);
        b = flagbit(obj,c);
        n = sum(bitset(0,b)); 
    end
    function c = code2flags(obj,n)
        validateattributes(n,{'numeric'},{'integer','scalar','nonnegative'});
        B = bitget(n,1:numel(obj.flags));
        c = obj.flags(B > 0); 
    end
   
    function C = merge(A,B,operation)
        if nargin < 3, operation = @bitor; end
        if numel(B) == 0, C = A; return; end
        if numel(A) == 0, C = B; return; end
                
        C = meteoQC();
        if isequal(A.flags,B.flags)
            placebitsfromA = @(X) X;
            placebitsfromB = @(X) X;
        else
            [ica,icb,ia,ib] = indexmembers(A.flags,B.flags);
            C.flags = [A.flags(ica),A.flags(ia),B.flags(ib)];
            placebitsfromA = @(X) reorderbits(X,[ica;ia],1:numel(A.flags));
            placebitsfromB = @(X) reorderbits(X,[icb;ib],[ica,numel(A.flags)+(1:numel(ib))]);
        end
        
        if isempty(A)
            [fica,ficb,fia,fib] = deal([],[],1:size(A,2),[]);
        elseif isempty(B)
            [fica,ficb,fia,fib] = deal([],[],[],1:size(B,2));
        else
            [fica,ficb,fia,fib] = indexmembers(fieldnames(A),fieldnames(B));
        end
        
        for j = 1:numel(fica)
            fld = A.data.Properties.VariableNames{fica(j)}; 
            valA = A.data{:,fica(j)};
            valB = B.data{:,ficb(j)};
            assert(isequal(size(valA),size(valB)),'Non matching field %s',fld);
            C.data.(fld) = operation(placebitsfromA(valA),placebitsfromB(valB));
        end
        for j = fia'
            fld = A.data.Properties.VariableNames{j};
            C.data.(fld) = operation(placebitsfromA(A.data{:,j}),0);
        end
        for j = fib'
            fld = B.data.Properties.VariableNames{j};
            C.data.(fld) = operation(placebitsfromB(B.data{:,j}),0);
        end
        
        function [ica,icb,ia,ib] = indexmembers(A,B)
        %   A(ica) = B(icb) are elemenents common to both A and B
        %   A(ia) are elements only found in A
        %   B(ib) are elements only found in B
        
            [~,ia,ib] = setxor(A,B);
            [~,ica,icb] = intersect(A,B);
        end
        
        function B = reorderbits(A,ia,ib)
        % make bitget(B,ib) = bitget(A,ia)
        
             % validateattributes(ia,{'numeric'},{'integer','positive'});
             % validateattributes(ib,{'numeric'},{'integer','positive','numel',numel(ia)});

             B = zeros(size(A),meteoQC.TYPE);
             for k = 1:numel(ia)
                 b = bitget(A,ia(k));
                 B = bitset(B,ib(k),b);
             end
        end
    end

    function C = mergeflags(F,old,new,operation)
        if nargin < 4, operation = @bitor; end
        
        old = cellstr(old);
        new = cellstr(new);
        validateattributes(new,{'cell'},{'scalar'},'','old');
        validateattributes(old,{'cell'},{'nonempty'},'','new');
        
        b = flagbit(F,old);
        old = F.flags(b);
        
        C = F;
        [ia,ib] = ismember(new,old);
        if ia == 0
            ib = 1;
            C.flags(b(1)) = new;
            % if isscalar(old), F = C; return; end
        end
        ir = [1:ib-1,ib+1:numel(b)];

        for j = 1:size(F,2)
           V = F.data{:,j};
           B = bitget(V,b(ib));
           for k = ir
               B = operation(B,bitget(V,b(k)));
           end
           C.data{:,j} = bitset(V,b(ib),B);
        end
        C = rmflag(C,b(ir));
    end
    
    function cat(d,varargin)
        
        validateattributes(d,{'numeric'},{'integer','scalar','positive','<',3});
        
        flaglist = cellfun(@(x) x.flags,varargin,'unif',0);
        flaglist = uniquecell([flaglist{:}],'stable');

        % bring to common flag encoding
        varargin = cellfun(@(x) setflags(x,flaglist),varargin,'unif',0);

        X = meteoQC;
        X.flags = flaglist;
        X.data = cat(d,varargin.data);
    end
    
    function A = setflags(A,flaglist)
        [ia,ib] = ismember(A.flags,flaglist);
        if isequal(ib,1:numel(ia)), return; end
        if ~all(ia)
            warning('Discarding information from %s',shortliststr(A.flags(~ia),'flag'))
        end
        
        X = A.data{:,:};
        B = zeros(size(X),meteoQC.TYPE);
        ib = ib(ia);
        ia = find(ia);
        for k = 1:numel(ia)
            B = bitset(B,ib(k),bitget(X,ia(k)));
        end
        A.data{:,:} = B;
        A.flags = flaglist;
    end
    
    function C = rmflag(F,b)
        C = F;
        n = numel(F.flags);
        if ~isnumeric(b), b = flagbit(F,b); end
        r = n-numel(b);
        for j = 1:size(C,2)
            V = C.data{:,j};
            N = numel(V);
            B = bitget(repmat(V(:),1,n),repmat(1:n,N,1));
            B(:,b) = [];
            C.data{:,j}(:) = sum(bitset(0,repmat(1:r,N,1),B,meteoQC.TYPE),2);
        end
        C.flags(b) = [];
    end
    
    function F = checkfield(F,MD,okfcn,bit,fields,flagifnan)
    % F = CHECKFIELD(F,MD,OKFCN,BIT,FLD) - set BIT of field(s) F.(FLD) of meteoQC object F if the
    %   corresponding field(s) of MD fail the given test, i.e. OKFCN( MD.(FLD) ) = 0
    % F = CHECKFIELD(F,MD,OKFCN,KEY,FLD)  - set bit FLAGBIT(F,KEY)
    % F = CHECKFIELD(F,MD,OKFCN,KEY,FLD,FLAGIFNAN)  - by default, NaN values in MD.(FLD) are not
    %   flagged. Use FLAGIFNAN = 1 when checking for missing values.
    %
    % See also: meteoQC.check, meteoQC.flagbit
    
        narginchk(4,6);
        if nargin < 5, fields = fieldnames(MD); end
        if nargin < 6, flagifnan = false; end
        
        if ~isnumeric(bit), [bit,F] = flagbit(F,bit); end
        validateattributes(bit,{'numeric'},{'integer','scalar','positive'});
        validateattributes(okfcn,{'function_handle'},{'scalar'});
        validateattributes(fields,{'char','cell','string'},{'nonempty'});
        
        if ~all(isfield(F,fields))
            [fields,idx] = parselist(fields,fieldnames(F),'-soft');
            fields = cellstr(fields);
            
            if any(idx == 0)
                % warning(shortliststr(fields(idx == 0),'Adding missing field','quotes',''''));
                for j = find(idx == 0)'
                    if isempty(j), continue; end
                    F.data.(fields{j}) = zeros(size(MD.(fields{j})),meteoQC.TYPE);
                end
            end
        else
            fields = cellstr(fields);
        end

        for j = 1:numel(fields)
           F.data.(fields{j}) = meteoQC.check(MD.(fields{j}),okfcn,bit,F.data.(fields{j}),flagifnan);
        end
    end

    function msg = flagsummary(S,FLD,filter,short)
    % Generate a summary string of the type:
    %
    %     'ENI: X/Y values flagged (0.0% num, 0.0% very_low, ...)
    %      GHI.1: X/Y values flagged (0.0% num, 0.0% weirdly_high, ...)
    %      ...'

        if nargin < 2 || (isempty(FLD) && ~iscell(FLD)), FLD = fieldnames(S); 
        else, FLD = parselist(FLD,fieldnames(S)); 
        end
        if nargin < 3, filter = true; end
        if nargin < 4, short = false; end
        
        Nflags = numel(S.flags);
        
        [NA_bit,~] = flagbit(S,'NA');
        
        FLD = cellstr(FLD);
        msg = repmat({{}},numel(FLD),1);
        for j = 1:numel(FLD)
            F = S.data.(FLD{j});
            [N,M] = size(F);
            if all(filter)
                Nf = N; 
                suffix = '';
            else
                Nf = nnz(filter);
                suffix = sprintf(' (from %0.1f%% subset)',Nf/N*100);
            end
            
            nbad = sum(F > 0,1) & filter;
            if all(nbad == 0)
                if ~short, msg{j} = {[FLD{j} ': no values flagged, 100% available' suffix]}; end
                continue; 
            end
            
            msg{j} = cell(M,1);
            for k = 1:M
    
                available = bitget(F(:,k),NA_bit) == 0 & filter;
                a = nnz(available);
                
                nbad = sum(F > 0,1) & filter;
            
                if M == 1, ID = FLD{j}; else, ID = sprintf('%s.%d',FLD{j},k); end
                if nbad(k) == 0
                    if ~short
                        msg{j}{k} = sprintf('%s: no values flagged, %0.1f%% available%s',...
                            ID,a/Nf*100,suffix); 
                    end
                    continue; 
                end
                msg{j}{k} = sprintf('%s: %d/%d values flagged, %0.1f%% available%s',...
                    ID,nbad(k),a,a/Nf*100,suffix);

                B = bitget(repmat(F(:,k),1,Nflags),repmat(1:Nflags,N,1));
                share = sum(B,1)/a*100;
                % if nnz(share) == 1
                %     msg{j} = sprintf('%s (%s)',msg{j},FLAGS{share > 0});
                % else
                % b = find(share > 0);
                [~,b] = maxk(share,4);
                shmsg = arrayfun(@(l) sprintf('%0.1f%% %s',share(l),S.flags{l}),b,'unif',0);
                msg{j}{k} = sprintf('%s (%s)',msg{j}{k},shortliststr(shmsg,'',5,'ellipsis','...'));
                % end
            end
        end
        msg = cat(1,msg{:});
        msg(cellfun(@isempty,msg)) = [];
        msg = strjoin(msg,newline());
    end
    
    function F = filterstructure(F,varargin)
    % bitwise FILTERSTRUCTURE (bits are set if result > 0)
        fld = fieldnames(F);
        r = [];
        for j = 1:numel(fld)
            X = F.data.(fld{j});
            n = ceil(log2(double(max(X,[],'all'))))+1;
            if n == -Inf
                if isempty(r)
                    r = numel(filterstructure(zeros(size(X,1),1),varargin{:}));
                end
                T.(fld{j}) = zeros([r,size(X,2)],meteoQC.TYPE);
            else
                B = single(bitget(repmat(X(:),1,n),repmat(1:n,numel(X),1)));
                B = reshape(B,[size(X),n]);
                b = filterstructure(B,varargin{:}) > 0;
                T.(fld{j}) = zeros(size(b,1),size(b,2),meteoQC.TYPE);
                for k = 1:n
                    T.(fld{j}) = bitset(T.(fld{j}),k,b(:,:,k),meteoQC.TYPE);
                end
            end
        end
        F.data = struct2table(T);
    end
end
end
