function [F,MD,msg] = irradianceQC(MD,t,sunel,varargin)

    % Definition of FLAGS and their order, to be consistent accross subfunctions
    %
    % BSRN: https://wiki.pangaea.de/wiki/BSRN_Toolbox#Quality_Check
    %
    %  	1 	Measurement falls below the physically possible limit
    %  	2 	Measurement exceeds the physically possible limit
    %  	4 	Measurement falls below the extremely rare limit
    %  	8 	Measurement exceeds the extremely rare limit
    %  	16 	Compared with a corresponding measurement the measurement is too low
    %  	32 	Compared with a corresponding measurement the measurement is too high 
    % 
    narginchk(2,Inf);
    
    TESTS = {'num','BSRN','CIE','clearsky','uncertainty'};
    FLAGS = {'num','abs_phys','rare_phys','flat','interp',...
         'BSRN_abs_lo','BSRN_abs_hi','BSRN_rare_lo','BSRN_rare_hi',...
         'CIE_lo','CIE_hi','CIE_out',...
         'CS_lo','CS_hi',...
         'UNC_hi','UNC_lo' };
    FLD = {'ENI','GHI','DHI','BNI','GTI','USW'};
    
    UNC = []; % placeholder
    
    flagbit({},FLAGS); % initialize
    
    [OPT,varargin] = getflagoptions(varargin,{'-independent'});
    OPT.sequence = {'num','BSRN','CIE'};
    OPT.valid = @(x) isfinite(x);
    OPT.flag2nan = {'num','BSRN_abs_lo','BSRN_abs_hi','abs_phys'};
    OPT.iterations = 1;
    OPT.initflags = []; % flags from previous IRRADIANCEQC calls
    
    % Parse field names, starting with defaults: OPT.GHI = 'GHI', OPT.DHI = 'DHI', ...
    fields = cell2struct(FLD',FLD');
    [fields,varargin,field_is_def] = getpairedoptions(varargin,fields);
    
    OPT = getpairedoptions(varargin,OPT,'restchk');

    OPT.sequence = parselist(OPT.sequence,TESTS,'test','sequence');
    OPT.flag2nan = parselist(OPT.flag2nan,FLAGS,'flag','flag2nan');
    if OPT.independent, OPT.iterations = 1; end

    % Arrange all fields into structure D, grouped by field types FLD
    [D,t] = getfields(MD,t,fields,field_is_def);
    N = numel(t);
    
    if nargin < 3, sunel = []; end
    if ~isempty(sunel)
        if ischar(sunel), sunel = MD.(sunel); end
        validateattributes(sunel,{'numeric'},{'real','vector','>=',-90,'<=',90,'size',[N,1]}); 
    end

    F = initflags(D,OPT.initflags);
    
    % Get/check/complete top of atmosphere irradiance 
    kr2 = solarposition.extraradiation(t,1); % (mean_sun_earth_dist/sun_earth_dist(t))^2
    if ~isfield(D,'ENI')
       D.ENI = kr2*1361.6; % see Kopp & Lean (2011) and solarposition.extraradiation.
    else
    % Allow both the old solar constant of 1366 W/m² as well as the lowest values from SORCE/TIM, 
    % with some margin for modeling error. Replace flagged values with own estimates.
       F.ENI = checkvalue(D.ENI,@isfinite,flagbit('num'),F.ENI,true);
       F.ENI = checkvalue(D.ENI,@(x) x < kr2*1369,flagbit('abs_phys'),F.ENI);
       F.ENI = checkvalue(D.ENI,@(x) x > kr2*1355,flagbit('abs_phys'),F.ENI);
       D.ENI(F.ENI > 0) = kr2(F.ENI > 0)*1361.6;
    end
    
    nanbits = flagbit(OPT.flag2nan);
    for k = 1:OPT.iterations
    for l = 1:numel(OPT.sequence)
        
        switch OPT.sequence{l}
            case 'num', F = checkfields(D,OPT.valid,flagbit('num'),F,[],true);
            case 'BSRN', F = BSRN_tests(D,sunel,F);
            case 'CIE', F = CIE_tests(D,sunel,F);
            case 'uncertainty', F = uncertainty_test(MD,sunel,UNC,F);
            case 'clearsky', F = clearsky_test(MD,sunel,UNC,F);
            otherwise
        end
        
        if ~OPT.independent, D = flagged2nan(D,F,nanbits); end
    end
    end
    
    if nargout > 1, MD = reassign(MD,D,fields); end
    
    if nargout > 2, msg = flagsummary(F,FLAGS); end
end
% 
% function [idx,list] = parselist(list,LIST,name)
% % Make sure list is a (case insensitive) subset of cellstr LIST
% % Read list = 'all' as LIST, and list = {'^',..} as setdiff(LIST,list(2:end))
% % Return position indices idx such that list = LIST(idx), and parsed list.
% 
%     if ischar(list) && strcmpi(list,'all')
%         list = 1:numel(LIST);
%     elseif ischar(list), list = {list};
%     end
% 
%     if iscellstr(list) || isstring(list)
%         inverted = strcmp(list{1},'^');
%         if inverted, list(1) = []; end
%         [ia,idx] = ismember(lower(list),lower(LIST));
%     elseif isnumeric(list)
%         idx = list;
%         ia = ismember(idx,1:numel(LIST));
%         inverted = false;
%     else
%         error('Expecting numeric, char, or cellstr %s',name);
%     end
%     assert(all(ia),'Unknown %s in %s',shortliststr(list(ia),'element'),name);
%     if inverted, idx = find(~ia); end
%     list = LIST(idx)';
% end

function [D,t] = getfields(MD,t,FLDS,field_is_default) 
% Parse struct/table/timetable MD, and for each field P in structure FLDS, arrange a copy:
%
%   D.(P) = [ MD.(p{1}),...,MD.(p{n}) ], where p{j} = FLDS.(P){j}
%
% In other words, FLDS is expected to be a dictionary that maps field names on MD to 
% 'standardized' field names fields(FLDS), horizontally concatenatenating several fields 
% whenever numel(FLDS.(P)) > 1.
%
% Fields FLDS.(P){j} not found in MD are simply not copied. A warning is issued only if the
% field was provided explicitly, i.e. ~FIELD_IS_DEFAULT.(P)

    P = fieldnames(FLDS);
    n = numel(P);
    
    switch class(MD)
        case 'timetable'
            if ~isempty(t), t0 = t; else, t0 = []; end
            [MD,t] = timetable2struct(MD);
            assert(isempty(t0) || isequal(t0,t),'Inconsistent time vector t and MD.t');
        case 'table'
            MD = table2struct(MD,'toscalar',true);
        case 'struct'
        otherwise
            error('Expecting structure or (time)table MD');
    end
    t = parsetime(t);
    N = numel(t);
    
    allflds = cell(1,n);
    for j = 1:n
        allflds{j} = FLDS.(P{j});
        if ~iscell(allflds{j}), allflds{j} = allflds(j); end
    end
    allflds = [allflds{:}];
    allflds(cellfun(@isempty,allflds)) = [];
    parsestruct(MD,[],'opt',allflds,'-n','-r','-m','size',[N,NaN]); 
    
    % Compile each type of sensor measurement into fields named according to FLD
    available = fieldnames(MD);

    MD = struct2cell(MD);

    for j = 1:numel(P)
        fld = FLDS.(P{j});
        [ia,ib] = ismember(fld,available);
        if ~field_is_default.(P{j}) && ~all(ia)
            warning('%s not found',shortliststr(fld(~ia),'field'));
        end
        if ~any(ia), continue; end
        D.(P{j}) = [MD{ib(ia)}];
    end
end

function MD = reassign(MD,D,FLDS) 
% Revert GETFIELDS, i.e. copy back (modified) values in D to their original fields in MD.

    P = fieldnames(D);
    available = fieldnames(MD);

    for j = 1:numel(P)
        fld = FLDS.(P{j});
        [ia,ib] = ismember(fld,available);
        switch nnz(ia)
        case 0
            MD.(P{j}) = D.(P{j});
            
        case 1
            fld = available{ib};
            MD.(fld) = D.(P{j});
            
        case size(D.(P{j}),2)
            fld = available(ib);
            r = 0;
            for k = 1:numel(ia)
                d = size(MD.(fld{k}),2);
                MD.(fld{k}) = D.(P{j})(:,r+1:r+d);
                r = r+d;
            end
            assert(r == size(D.(P{j}),2),'Unexpected reassign condition 1');
        otherwise
            error('Unexpected reassign condition 2');
        end
    end
end

function F = BSRN_tests(MD,sunel,F)
% BSRN recoomended QC tests (absolute upper and lower bounds)
% https://bsrn.awi.de/fileadmin/user_upload/bsrn.awi.de/Publications/BSRN_recommended_QC_tests_V2.pdf

    if nargin < 2 || isempty(sunel), cz = 1; 
    else, cz = max(0,sind(sunel));
    end

    B = flagbit({'BSRN_abs_lo','BSRN_abs_hi','BSRN_rare_lo','BSRN_rare_hi'});

    F = checkfields(MD,@(x) x > -4.0,B(1),F);
    F = checkfields(MD,@(x) x > -2.0,B(3),F);

    if isfield(MD,'GHI')  
        F.GHI = checkvalue(MD.GHI,@(x) x < 1.5*MD.ENI.*cz.^1.2 + 100,B(2),F.GHI);
        F.GHI = checkvalue(MD.GHI,@(x) x < 1.2*MD.ENI.*cz.^1.2 + 50,B(4),F.GHI);
    end

    if isfield(MD,'DHI')
        F.DHI = checkvalue(MD.DHI,@(x) x < 0.95*MD.ENI.*cz.^1.2 + 50,B(2),F.DHI);
        F.DHI = checkvalue(MD.DHI,@(x) x < 0.75*MD.ENI.*cz.^1.2 + 30,B(4),F.DHI);
    end

    if isfield(MD,'BNI')      
        F.BNI = checkvalue(MD.BNI,@(x) x < MD.ENI,B(2),F.BNI);
        F.BNI = checkvalue(MD.BNI,@(x) x < 0.95*MD.ENI.*cz.^0.2 + 10,B(4),F.BNI);
    end
    
    if isfield(MD,'USW')
        F.USW = checkvalue(MD.USW,@(x) x < 1.2*MD.ENI.*cz.^1.2 + 50,B(2),F.USW);
        F.USW = checkvalue(MD.USW,@(x) x < MD.ENI.*cz.^1.2 + 50,B(4),F.USW);
    end
    
    % There is no explicit recommendation for GTI, so this is just GHI + BNI
    % TODO: consider incidence angle
    if isfield(MD,'GTI')
        F.GTI = checkvalue(MD.GTI,@(x) x < MD.ENI.*(0.95*cz.^1.2 + 1) + 50,B(2),F.GTI);
        F.GTI = checkvalue(MD.GTI,@(x) x < MD.ENI.*(0.75*cz + 0.95).*cz.^0.2 + 40,B(4),F.GTI);
    end
end

function F = CIE_tests(MD,sunel,F)
% CIE recommended QC tests: closure within 15%, DHI < GHI within 10% flagged as CIE_hi / CIE_lo
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

    if ~isfield(MD,'GHI'), return; end

    if nargin < 2 || isempty(sunel)
        cz = [];
        sunel = [];
    else
        cz = max(0,sind(sunel));
    end

    B = flagbit({'CIE_lo','CIE_hi','CIE_out'});
    
    G0 = max(0,max(MD.GHI,[],2));
    if ~isempty(cz) && all(isfield(MD,{'BNI','DHI'}))
    % Get all combinations of BHI + DHI, arrange in dimensions {Nt,GHI,BNI,DHI}
        GHIc = permute(MD.BNI.*cz,[1,3,2]) + permute(MD.DHI,[1,3,4,2]);
        G0 = max(G0,max(GHIc,[],3:4));
    else
        GHIc = [];
    end
    
    if isempty(sunel)
        inrange = G0 >= 20;
        warning('Possibly applying CIE QC checks outside recommended conditions'); 
    else
    % Flag values only if sunel > 4 and GHI > 20 W/m²
        inrange = sunel >= 4 & G0 >= 20;
    end
    F = checkfields(MD,@(x) inrange,B(3),F);

    if isfield(MD,'DHI')
    % Check GHI - DHI > 10% GHI, flag GHI and DHI as lo/hi otherwise

        LO_global = MD.GHI < permute(MD.DHI,[1,3,2]) - G0*0.1;
        notlo = ~all(LO_global,3);
        F.GHI = checkvalue(MD.GHI,@(x) ~inrange | notlo,B(1),F.GHI);

        nothi = permute(~all(LO_global,2),[1,3,2]);
        F.DHI = checkvalue(MD.DHI,@(x) ~inrange | nothi,B(2),F.DHI);
    end

    if ~isempty(GHIc)
    % Check closure: | BNI cz + DHI - GHI | < 15% GHI, flag each as lo/hi otherwise
    
        na = isnan(MD.GHI) | isnan(GHIc);
        LO_global = MD.GHI < 0.85*GHIc | na; % {Nt,GHI,BNI,DHI}
        HI_global = MD.GHI > 1.15*GHIc | na;
        
        notlo = ~all(LO_global,3:4);
        nothi = ~all(HI_global,3:4);
        F.GHI = checkvalue(MD.GHI,@(x) ~inrange | notlo,B(1),F.GHI);
        F.GHI = checkvalue(MD.GHI,@(x) ~inrange | nothi,B(2),F.GHI);

        notlo = permute(~all(HI_global,2:3),[1,4,2,3]);
        nothi = permute(~all(LO_global,2:3),[1,4,2,3]);
        F.DHI = checkvalue(MD.DHI,@(x) ~inrange | notlo,B(1),F.DHI);
        F.DHI = checkvalue(MD.DHI,@(x) ~inrange | nothi,B(2),F.DHI);

        notlo = permute(~all(HI_global,[2,4]),[1,3,2,4]);
        nothi = permute(~all(LO_global,[2,4]),[1,3,2,4]);
        F.BNI = checkvalue(MD.BNI,@(x) ~inrange | notlo,B(1),F.BNI);
        F.BNI = checkvalue(MD.BNI,@(x) ~inrange | nothi,B(2),F.BNI);
    end
end
    
function F = uncertainty_test(MD,sunel,UNC,F)
% EXPERIMENTAL: Uncertainty-based QC

    B = flagbit({'UNC_lo','UNC_hi'});
    
    [MD0,UNC] = combinedguess(MD,UNC,sind(sunel));
    
    fields = intersect(fieldnames(MD),fieldnames(UNC));
    for j = 1:numel(fields)
        notlo = MD.(fields{j}) >= MD0.(fields{j}) - UNC.(fields{j});
        F.(fields{j}) = checkvalue(MD.(fields{j}),@(x) notlo,B(1),F.(fields{j}));
        
        nothi = MD.(fields{j}) <= MD0.(fields{j}) + UNC.(fields{j});
        F.(fields{j}) = checkvalue(MD.(fields{j}),@(x) nothi,B(2),F.(fields{j}));
    end
end

function F = clearsky_test(MD,sunel,UNC,F)
% EXPERIMENTAL: Clear-Sky-model-based QC

    REQ = {'CSGHI','CSBNI','CSDHI','ENI','CS'};
    DEBUG = true;
    MINCSFRACTION = 0.05;
    B = flagbit({'CS_lo','CS_hi'});
    
    missing = ~isfield(MD,REQ);
    if any(missing)
        warning('Skipping clearsky test. %s',shortliststr(REQ(missing),'Missing field'));
        return;
    end
    
    dark = sunel < 0 & MD.CSGHI > 0;
    if all(isfield(MD,{'DHI','GHI'})), dark = dark | ~(MD.GHI > 0) | ~(MD.DHI > 0); 
    elseif isfield(MD,'BNI'), dark = dark | ~(MD.BNI > 0);
    else, return;
    end
    
    if nnz(MD.CS(~dark))/nnz(~dark) < MINCSFRACTION
        warning('Not enough clear samples for clear-sky test');
        return;
    end

    % Use best available (inverse-variance-weighted) estimate for BNI up to now.
    [CG,CU] = combinedguess(MD,UNC,sind(sunel));

    if isfield(MD,'BNI')
    % Check that BNI < CSBNI
        bad = MD.BNI > MD.CSBNI;
        if any(bad,'all')
            
            % Get an estimate of the minimum CSkn - kn, based on "sunny" samples ...
            % "sunny" =  clear-sky | any point with BNI > CSBNI
            f = MD.CS | CG.BNI > MD.CSBNI & ~dark;  
            CSkn = MD.CSBNI./MD.ENI;
            kn = CG.BNI(f)./MD.ENI(f);
            sn = min(CU.BNI(f)./MD.ENI(f));
            [lb,~,out] = binnedIQR(CSkn(f),(CSkn(f)-kn),'-asym','minrange',sn,'xlim',[0 1],'-interp');
            kn_max = CSkn - lb(CSkn);
            bad = BNI./MD.ENI > kn_max; % ... to use as min. threshold
            
            F.BNI = checkvalue(MD.BNI,@(x) ~bad,B(2),F.BNI);
            MD.CS = MD.CS & ~out;
            
            if DEBUG
                GUIfigure('CSBNI','Clear-Sky BNI check'); clf(); hold on; % #ok<UNRCH>
                CSkn = repmat(CSkn,1,size(BNI,2));
                scatter(CSkn(:),reshape(BNI./MD.ENI,[],1),1,[1 1 1]*0.8);
                densityplot(CSkn(f,1),kn(~out),2);
                plot(lb.GridVectors{1},lb.GridVectors{1}-lb.Values,'r:')
                kn = BNI./MD.ENI;
                scatter(CSkn(bad),kn(bad),2,'r');      
                grid on; plot([0,1],[0,1],'k:'); 
                axis equal; axis([0 0.8 0 0.8]);
                xlabel('Clear-Sky kn = BNI_{CS}/G_{extra}');
                ylabel('Measured kn = BNI_{CS}/G_{extra}');
                clear kn CSkn sn out
            end
        end
    end

    if all(isfield(MD,{'DHI','GHI'}))
        % Check that CSkd < kd
        CSkd = MD.CSDHI./MD.CSGHI;
        bad = any(MD.DHI < CSkd.*CG.GHI,2) | any(CG.DHI < MD.GHI.*CSkd,2);
        if any(bad)
        % Get the minimum fraction kd/CSkd based on clear-sky samples
    
            f = MD.CS;
            kd = CG.DHI(f)./CG.GHI(f);
            sd = sqrt((kd.*CU.GHI(f)).^2 + CU.DHI(f).^2)./CG.GHI(f); % conservative: assumes no correlation
            [lb,~,out] = binnedIQR(CSkd(f),kd./sd./CSkd(f),'-asym','xlim',[0 1],'minrange',1,'-interp');
            f(f) = ~out;
            kd =[MD.DHI./CG.GHI,CG.DHI./MD.GHI];
            kd_min = lb(CSkd).*CSkd; % ... minimum expected CSkd

            F.DHI = checkvalue(MD.DHI,@(x) MD.DHI./CG.GHI > kd_min,B(1),F.DHI);
            F.GHI = checkvalue(MD.GHI,@(x) CG.DHI./MD.GHI > kd_min,B(2),F.GHI);

            if DEBUG
                GUIfigure('CSkd','Clear-Sky kd check'); clf(); hold on; % #ok<UNRCH>
                CSkd = repmat(CSkd,1,size(kd,2));
                scatter(CSkd(:),kd(:),1,[1 1 1]*0.8);
                densityplot(CSkd(f,1),CG.DHI(f)./CG.GHI(f),2);
                plot(lb.GridVectors{1},lb.Values.*lb.GridVectors{1},'r:')
                scatter(CSkd(bad),kd(bad),2,'r');
                grid on; plot([0,1],[0,1],'k:');
                axis equal; axis([0 0.8 0 0.8]);
                xlabel('Clear-Sky kd = DHI_{CS}/GHI_{CS}');
                ylabel('Measured kd = DHI/GHI*, DHI*/GHI')
                clear kd sd out CSkd
            end
        end
    end
end

function b = flagbit(flag_keys,flag_list)
% flagbit(~,list) - initialize: set persistent flag list
% flagbit(keys) - return indices in in flag list corresponding to keys

    persistent FLAGS;
    if nargin > 1, FLAGS = flag_list; end
    if isempty(FLAGS), FLAGS = flag_keys; end % initialize on first call
    
    if nargin == 0 || (~iscell(flag_keys) && isempty(flag_keys))
        b = 1:numel(FLAGS);
        return;
    end
    if isempty(flag), b = []; return; end
    if ~iscell(flag_keys), flag_keys = {flag_keys}; end

    [ia,b] = ismember(flag_keys,FLAGS);
    assert(all(ia),'Unknown flag key: %s',flag_keys(ia));
end

function F = initflags(D,F0)
% Create a structure similar to D, with uint32 flag values (set to 0)

    fld = fieldnames(D);

    if isempty(F0), F0 = struct(); else
        F0 = parsestruct(F0,{},'opt',fld,'-n','-i','-p','-reduced');
    end

    F = struct();
    for j = 1:numel(fld)
        if isfield(F0,fld)
            assert(isequal(size(F0.(fld{j})),size(D.(fld{j}))),'Incompatible start flags');
            F.(fld{j}) = F0.(fld{j});
        else
            F.(fld{j}) = zeros(size(D.(fld{j})),'uint32');
        end
    end
end

function flag = checkvalue(value,okfcn,key,flag,flagifnan)
% F = check(V,fcn,B,F) - if ~fcn(V) & isfinite(V), set bit B on F
% F = check(V,fcn,KEY,F) - ... set bit corresponding to KEY on F
% F = check(V,fcn,B,F,1) - set B  if ~fcn(V), even if ~isfinite(V)

    if nargin < 5, flagifnan = false; end
    if ~isnumeric(key), key = MeteoData.flagbit(key); end

    notok = ~okfcn(value);
    if ~flagifnan
        notok = notok & ~isnan(value);
    end
    flag(notok) = bitset(flag(notok),key);
end 

function F = checkfields(MD,okfcn,key,F,fields,flagifnan)
% Do F.(f) = check(S.(f),fcn,key,F.(f)) for each FIELD

    if nargin < 5 || (~iscell(fields) && isempty(fields)), fields = fieldnames(F); end
    if nargin < 6, flagifnan = false; end

   for j = 1:numel(fields)
       F.(fields{j}) = checkvalue(MD.(fields{j}),okfcn,key,F.(fields{j}),flagifnan);
   end
end

function X = flagged2nan(X,F,flag_bits,fields)
% B = flagged2nan(A,F,flag_bits) - replace elements in A for which any(bitget(F,flag_bits)) by NaN
% D = flagged2nan(D,F,flag_bits,fld) - do the same for every field fld of structure D.

    if isstruct(X) && (nargin < 4 || (~iscell(fields) && isempty(fields)))
        fields = fieldnames(F);
    end

    for k = 1:numel(flag_bits)
        if isstruct(X)
            for j = 1:numel(fields)
                notok = logical(bitget(F.(fields{j}),flag_bits(k)));
                X.(fields{j})(notok) = NaN;
            end
        else
            notok = logical(bitget(F,flag_bits(k)));
            X(notok) = NaN;
        end
    end
end

function msg = flagsummary(S,FLAGS)
% Generate a summary string of the type:
%
%     'ENI: X/Y values flagged (0.0% num, 0.0% rare_lo, ...)
%      GHI: X/Y values flagged (0.0% num, 0.0% rare_lo, ...)
%      ...'

    FLD = fieldnames(S);
    msg = cell(numel(FLD),1);
    for j = 1:numel(FLD)
        F = S.(FLD{j});
        [N,M] = size(F);

        Nf = nnz(F > 0);
        if Nf == 0, msg{j} = [FLD{j} ': no values flagged']; continue; end

        msg{j} = sprintf('%s: %d/%d values flagged',FLD{j},Nf,N*M);
        if size(F,2) > 1
            r = any(F > 0,1);
            msg{j} = sprintf('%s on %d/%d channels',msg{j},nnz(r),M);
        end

        % Expand F bit-wise, so that F(i,j,k) is the jth bit of observation i, variable k
        B = fliplr(dec2bin(F(:),numel(FLAGS)) == '1'); % LSB to MSB
        B = permute(reshape(B,N,M,[]),[1,3,2]);

        share = sum(B,[1,3])/(N*M);
%         if nnz(share) == 1
%             msg{j} = sprintf('%s (%s)',msg{j},FLAGS{share > 0});
%         else
            share = arrayfun(@(j) sprintf('%0.1f%% %s',share(j)*100,FLAGS{j}),find(share > 0),'unif',0);
            msg{j} = sprintf('%s (%s)',msg{j},strjoin(share,', '));
%         end
    end
    msg = strjoin(msg,newline());
end

% %% Sensor merge/closure, kt-kd check, and gap-filling
%     
%     % if (noBNI && noDHI), break; end
%     
%     % Update best estimates (joint MLE of all sensors)
%     [g0,d0,b0] = combinedguess(GHI,DHI,BNI,UNC,cosZ);
% 
%     % Look for points that do not match best estimates with probability > opt.P
%     K = abs(norminv(opt.P))/2; % assuming UNC contains 95% 'expanded' uncertainties
% 
%     bad = abs(GHI-g0) > UNC.GHI.*K;
%     GHI(bad) = NaN;
%     meh = any(bad,2);
% 
%     bad = abs(DHI-d0) > UNC.DHI.*K;
%     DHI(bad) = NaN;
%     meh = meh | any(bad,2);
% 
%     bad = abs(BNI-b0) > UNC.BNI.*K;
%     DHI(bad) = NaN;
%     meh = meh | any(bad,2);
% 
%     if any(meh) && QCpass == 1   
%         fixreport('Cross-check',:) = {nnz(meh & available),nmbe(BNI,b0),nmbe(DHI,d0),nmbe(GHI,g0)};
%         last.GHI = GHI;
%         last.DHI = DHI;
%         last.BNI = BNI;
%     end 
% 
%     % Update best joint estimates (for good)
%     [GHI,DHI,BNI,UNC.GHI,UNC.DHI,UNC.BNI] = combinedguess(GHI,DHI,BNI,UNC,cosZ);
% 
%     if size(last.GHI,2) > 1, last.GHI = mean(last.GHI,2,'omitnan'); end
%     if size(last.DHI,2) > 1, last.DHI = mean(last.DHI,2,'omitnan'); end
%     if size(last.BNI,2) > 1, last.BNI = mean(last.BNI,2,'omitnan'); end
% 
%     meh = any(abs([last.GHI,last.DHI,last.BNI]-[GHI,DHI,BNI])./max([GHI,DHI,BNI],opt.minGHI) > opt.reltol,2);
%     if any(meh) && QCpass == 1
%         fixreport('Closure',:) = {nnz(meh & available),nmbe(BNI,last.BNI),nmbe(DHI,last.DHI),nmbe(GHI,last.GHI)};
%         last.GHI = GHI;
%         last.DHI = DHI;
%         last.BNI = BNI;
%     end
% 
%     % Find abnormal kd-kt pairs...
%     % TODO: merge into (soft) max. likelihood approach, applying a penalty based on kd-kt plot
%     %   kernel density for "normal" data.
%     kd = DHI./GHI;
%     kt = GHI./EHI; 
%     kt(dark) = NaN;
%     kd = min(1,kd);
%     kd(dark) = NaN;
%     kd(~(GHI > 0)) = NaN;
%     bad = ktkdfilter(kt,kd,'-silent');
%     bad = bad & isfinite(kt) & isfinite(kd);
% 
%     if any(bad,'all')
%         GHI(bad) = NaN; DHI(bad) = NaN; BNI(bad) = NaN; kt(bad) = NaN; kd(bad) = NaN;
%         if QCpass == 1
%             fixreport('kt-kd',:) = {nnz(bad & available),nmbe(BNI,last.BNI),nmbe(DHI,last.DHI),nmbe(GHI,last.GHI)};
%             last.GHI = GHI; last.DHI = DHI; last.BNI = BNI;
%         end
%     end
% 
%     if QCpass > 1 && nmbe(BNI,last.BNI) == 0 && nmbe(DHI,last.DHI) == 0 && nmbe(GHI,last.GHI) == 0
%         break;
%     end
% end
% 
%     meh = any(abs([last.GHI,last.DHI,last.BNI]-[GHI,DHI,BNI])./max([GHI,DHI,BNI],opt.minGHI) > opt.reltol,2);
%     if any(meh)
%         fixreport('Re-checks',:) = {nnz(meh & available),nmbe(BNI,last.BNI),nmbe(DHI,last.DHI),nmbe(GHI,last.GHI)};
%         last.GHI = GHI;
%         last.DHI = DHI;
%         last.BNI = BNI;
%     end
% 
%     tweaked.GHI = abs(MD.GHI - GHI) > max(eps(GHI),opt.reltol.*GHI);
%     tweaked.DHI = abs(MD.DHI - DHI) > max(eps(DHI),opt.reltol.*DHI);
%     tweaked.BNI = abs(MD.BNI - BNI) > max(eps(BNI),opt.reltol.*BNI);
%     flagged.GHI = ~isfinite(GHI) & isfinite(MD.GHI);
%     flagged.DHI = ~isfinite(DHI) & isfinite(MD.DHI);
%     flagged.BNI = ~isfinite(BNI) & isfinite(MD.BNI);
%     
%     tofill = ~(isfinite(MD.GHI) & isfinite(MD.DHI) & isfinite(MD.BNI) & ~dark);
%     
%     if opt.fillgaps > 0 && any(tofill)
%     % Interpolate NaN values in irradiance
%  
%         % Fill gaps using coupled clear-sky indices
%         Kc = [GHI./MD.CSGHI,DHI./MD.CSDHI,BNI./MD.CSBNI];
%         Kc(dark,:) = NaN;
%         [~,Kc] = guessmissing(Kc,dt,0.25,'',missing);
%         
%         % Let any gap wider than 2·opt.fillgaps as NaN
%         semiavailable = isfinite(MD.GHI) + isfinite(MD.DHI) + isfinite(MD.BNI) >= 2;
%         idx = find(semiavailable);
%         d = interp1(idx,idx,1:Nt,'nearest','extrap')-(1:Nt); % distance to nearest valid point
%         Kc(abs(d) > opt.fillgaps,:) = NaN;
%         
%         g0 = Kc(:,1).*MD.CSGHI; g0(dark) = 0;
%         d0 = Kc(:,2).*MD.CSDHI; d0(dark) = 0;
%         b0 = Kc(:,3).*MD.CSBNI; b0(dark) = 0;
%         U.GHI = rms(g0(~dark) - GHI(~dark));
%         U.DHI = rms(d0(~dark) - DHI(~dark));
%         U.BNI = rms(b0(~dark) - BNI(~dark));
%         [g0,d0,b0] = combinedguess(g0,d0,b0,U,cosZ);
%         
%         bad = ktkdfilter(g0./EHI,d0./g0,'-silent');
%         g0(bad) = NaN; d0(bad) = NaN; b0(bad) = NaN;
% 
%         f = ~isfinite(GHI) & isfinite(g0);
%         if any(f)
%             GHI(f) = g0(f);
%             UNC.GHI(f) = Inf;
%             warning('cmd:missing','Interpolating %d/%d missing points in GHI - REALLY NOT RECOMMENDED!',nnz(f),numel(f));
%         end
%         f = ~isfinite(DHI) & isfinite(d0);
%         if any(f)
%             DHI(f) = d0(f);
%             UNC.DHI(f) = Inf;
%             warning('cmd:missing','Interpolating %d/%d missing points in DHI - REALLY NOT RECOMMENDED!',nnz(f),numel(f));
%         end
%         f = ~isfinite(BNI) & isfinite(b0);
%         if any(f)
%             BNI(f) = b0(f);
%             UNC.BNI(f) = Inf;
%             warning('cmd:missing','Interpolating %d/%d missing points in BNI - REALLY NOT RECOMMENDED!',nnz(f),numel(f));
%         end
%         
%         fixreport('Gap-filling',:) = {nnz(tofill),nmbe(BNI,last.BNI),nmbe(DHI,last.DHI),nmbe(GHI,last.GHI)};
% 
%         kd = DHI./GHI;
%         kt = GHI./EHI; 
%         kt(dark) = NaN;
%         kd = min(1,kd);
%         kd(dark) = NaN;
%         kd(~(GHI > 0)) = NaN;
%         
%     elseif strcmpi(opt.withintol,'nan')
%         GHI(tweaked.GHI) = NaN; 
%         DHI(tweaked.DHI) = NaN; 
%         BNI(tweaked.BNI) = NaN;
%         bad = tweaked.GHI | tweaked.BNI | tweaked.DHI;
%         fixreport('tweaked = NaN',:) = {nnz(bad),nmbe(BNI,last.BNI),nmbe(DHI,last.DHI),nmbe(GHI,last.GHI)};
%     end
%     
%     f = [nmbe(MD.BNI,BNI),nmbe(MD.DHI,DHI),nmbe(MD.GHI,GHI)];
%     f(2,:) = [nrms(MD.BNI,BNI),nrms(MD.DHI,DHI),nrms(MD.GHI,GHI)];
%     if ~isempty(fixreport) || any(abs(f),'all') > opt.reltol
%         n = nnz((tweaked.GHI | tweaked.BNI | tweaked.DHI | flagged.GHI | flagged.BNI | flagged.DHI) & available);
%         fixreport('total (nMBE)',:) = {n,f(1,1),f(1,2),f(1,3)};
%         fixreport('total (nRMS)',:) = {n,f(2,1),f(2,2),f(2,3)};
%     end
% 
%     % Put everything back into place
%     % NaNs have been placed (if required) in GHI,.. already
%     MD.tweaked = tweaked;
%     MD.flagged = flagged;
%     MD.uncertainty = UNC;
%     if any(strcmpi(opt.withintol,{'fix','nan'}))  
%         MD.GHI = GHI;
%         MD.BNI = BNI;
%         MD.DHI = DHI;
%         MD.tweaked = tweaked;
%         printif('\nThe following changes have been applied after Quality-Control\n');
%     else
%         warning('Quality-Control adjustments have been detected, but NOT APPLIED!');
%     end
%     if any(strcmpi(opt.beyondtol,{'fix','nan'}))
%         MD.GHI(flagged.GHI) = NaN;
%         MD.DHI(flagged.DHI) = NaN;
%         MD.BNI(flagged.BNI) = NaN;
%     end
%     MD.kt = kt; % best estimates for kd,kt stay, no matter what
%     MD.kd = kd;
% 
% %% User feedback
% 
%     if ~isempty(fixreport)
%         fixreport{:,'F'} = fixreport{:,'F'}./nnz(available);
%         fixreport{:,:} = round(fixreport{:,:}*100,2);
% 
%         printif('[ correction type, fraction of samples afected, effect on BNI, DHI, GHI ]\n');
%         printif('(unless stated, values are MBE percentages based on available daylight samples):\n\n');
%         if opt.verbose, disp(fixreport); end
%     end
%     
%     meh = false; bad = false;
% 
%     f = nnz((flagged.BNI | flagged.DHI | flagged.GHI) & available)/nnz(available);
%     if f > opt.reltol
%         msg = {sprintf('%0.2f%% of daylight points flagged as incorrect/outliers',f*100)};
%         meh = f > opt.warning;
%         bad = f > opt.error;
%         if ~(bad || meh), printif('%s\n',msg{end}); msg(end) = []; end % print and dump
%     else
%         msg = {};
%     end
%     
%     f = nnz((tweaked.GHI | tweaked.BNI | tweaked.DHI) & available)/nnz(available);
%     if f > opt.reltol
%     	msg{end+1} = sprintf(['%0.2f%% of daylight points adjusted (within tolerance) ',...
%                 'to comply with physical limits & closure conditions'],f*100);
%         meh = meh || f > opt.warning;
%         % bad = bad || f > opt.error;
%         if ~(bad || meh), printif('%s\n',msg{end}); msg(end) = []; end
%     end
% 
%     f = fixreport{'total (nRMS)',2:end}/100;
%     if any(f) > opt.reltol 
%         if ~isempty(msg)
%             msg = sprintf('%s, resulting in changes (nRMS) of ',shortliststr(msg));
%         else
%             msg = 'Quality-Control checks result in changes (nRMS) of ';
%         end
%         msg = sprintf('%s%0.2f%% BNI, %0.2f%% DHI, and %0.2f%% GHI. ',msg,f*100);
%     else
%         if ~isempty(msg)
%             msg = sprintf('%s, bias falls within %0.1e tolerance. ',shortliststr(msg),opt.reltol);
%         else
%             msg = sprintf('Quality-Control checks fall within %0.1e tolerance. ',opt.reltol);
%         end
%     end
%     meh = meh || any(f > min(opt.warning,opt.tolerance));
%     bad = bad || any(f > opt.error);
%     if ~(bad || meh), printif('%s\n',msg); msg = ''; end
% 
%     available = isfinite(MD.GHI) & isfinite(MD.DHI) & isfinite(MD.BNI) & ~dark;
%     f = nnz(available)/nnz(~dark);
%     msg = sprintf('%sData availability within provided period: %0.2f%%.',msg,f*100);
%     meh = meh | 1-f > opt.warning;
%     if ~(bad || meh), printif('%s\n',msg); msg = ''; end
%     
%     % Issue warning/error based on QC status
%     if bad || meh, warning('cmd:QC','%s',msg); end
%     printif('\n');
% 
%     MD.available = available;  

function [x0,u0] = bestguess(x,u,dim)
% Inverse-variance weighted X, according to uncertainty U

    if nargin < 3, dim = 2; end
    [x,u] = compatiblesize(x,u);
    if size(x,dim) == 1, x0 = x; u0 = u; return; end
    
    w = 1./u.^2;
    w(isnan(x)) = NaN;
    u0 = 1./sum(w,dim,'omitnan');
    x0 = sum(x.*w,dim,'omitnan').*u0;
    u0 = sqrt(u0);
end

function [CG,CU] = combinedguess(MD,UNC,cosZ)
% Given a set of measurements GHI, DHI, BNI with uncertainties UNC.GHI, UNC.DHI, UNC.BNI; get a 
% maximum-likelihood estimate of the 'true' G0, D0, B0. The estimates must satisfy the closure 
% condition GHI = DHI + BNI·cosZ while maximizing the log-likelyhood:
%
%       L = ln{P(G0-GHI)·P(B0 - BNI)·P(D0-DHI)}
%
% FIX: Replace by a flexible State-Space Recursive Bayesian Estimator:
%   + Consider tilted sensors.
%   + Account for correlation between samples (sensor calibration & offset), and between sensors
%   (e.g. cosine-response, spectral-response, etc.), 
%   + Consider expected kt-kd distribution, e.g. from probabilistic separation model
%   + Account for uncertainty in closure equation for non-finite (downsampled) intervals!

    % FIX: currently not handling anything other than GHI, BNI, DHI
    CG = MD;
    CU = UNC; 

    % Start with inverse-variance weighted estimates
    if isfield(MD,'GHI'), [CG.GHI,CU.GHI] = bestguess(MD.GHI,UNC.GHI); end
    if isfield(MD,'DHI'), [CG.DHI,CU.DHI] = bestguess(MD.DHI,UNC.DHI); end
    if isfield(MD,'BNI'), [CG.BNI,CU.BNI] = bestguess(MD.BNI,UNC.BNI); end

    if all(isfield(MD,{'GHI','DHI','BNI'}))
    % For the simple case of uncorrelated variables the solution to the constrained MLE problem
    % can be calculated analytically:
    
        Ug = CU.GHI.^2; Ud = CU.DHI.^2; Ub = CU.BNI.^2;
        lambda = (CG.GHI - CG.DHI - cosZ.*CG.BNI)./(Ud + Ub.*cosZ.^2 + Ug);
        lambda(~isfinite(lambda)) = 0;

        CG.BNI = CG.BNI + lambda.*Ub.*cosZ;
    
        % Check and enforce constraint B.BNI > 0
        bad = CG.BNI < 0;
        lambda(bad) = (CG.GHI(bad) - CG.DHI(bad))./(Ud(bad) + Ug(bad));
        CG.BNI(bad) = 0;
    
        CG.GHI = CG.GHI - lambda.*Ug;
        CG.DHI = CG.DHI + lambda.*Ud;
        
    elseif all(isfield(MD,{'GHI','DHI'}))
        
        overcast = MD.GHI < MD.DHI;
        [x0,u0] = bestguess([CG.GHI(overcast),CG.DHI(overcast)],[CU.GHI(overcast),CU.DHI(overcast)]);
        CG.GHI(overcast) = x0;
        CG.DHI(overcast) = x0;
        CU.GHI(overcast) = u0;
        CU.DHI(overcast) = u0;
        
    end
   
    % % Check for closure: BNI·cos(Z) + DHI = GHI, or equivalently: kt(1-kd) = kn
    % % meh = abs(B0.*cosZ + D0 - G0) > opt.reltol.*EHI;
    % meh = all(isfinite([GHI,DHI,BNI,UNC.GHI,UNC.DHI,UNC.BNI]),2);
    % 
    % if any(meh)
    %     bad = meh;
    % 
    %     % Divide into chunks so that LSQNONLIN can work with reasonably-sized matrices
    %     N = ceil(nnz(meh)/round(nnz(meh)/200));  % ~200 pts seems to be close to the sweet spot
    %     i = 0; n = 0;
    %     for j = 1:ceil(nnz(meh)/N)
    %         i = find(meh(i(end)+1:end),N)+i(end);
    %         gi = GHI(i,:); bi = BNI(i,:); di = DHI(i,:);
    %         sg = UNC.GHI(i,:); sd = UNC.DHI(i,:); sb = UNC.BNI(i,:);
    %         cz = cosZ(i);
    % 
    %         % Log-likelihood function
    %         L = @(BD) [(BD(:,1).*cz+BD(:,2)-gi)./sg,(BD(:,2)-di)./sd,(BD(:,1)-bi)./sb];
    %         x0 = [B0(i),D0(i)];
    %         op = optimoptions('lsqnonlin','StepTolerance',opt.reltol,'Display','off');
    %         if n == numel(i)
    %         % J = full(sparse(mod(0:4*n-1,3*n)+1,mod(0:4*n-1,2*n)+1,1));
    %             op.JacobPattern = (J > 0)*1;
    %         end
    %         n = numel(i);
    %         [x0,~,y0,flag,~,~,J] = lsqnonlin(L,x0,[],[],op);
    % 
    %         bad(i) = prod(normcdf(-abs(y0))*2,2) < opt.P; % best joint likelyhood
    % 
    %         B0(i) = x0(:,1);
    %         D0(i) = x0(:,2);
    %         G0(i) = B0(i).*cz + D0(i);
    % 
    %         p = @(b,d) -sum(L([b,d]).^2,2)/2;
    % 
    %         db = 1; % Sb(i)/2;
    %         dd = 1; % Sd(i)/2;
    % 
    %         H = [p(B0(i)-db,D0(i)-dd),p(B0(i),D0(i)-dd),p(B0(i)+db,D0(i)-dd),...
    %              p(B0(i)-db,D0(i)   ),p(B0(i),D0(i)   ),p(B0(i)+db,D0(i)   ),...
    %              p(B0(i)-db,D0(i)+dd),p(B0(i),D0(i)+dd),p(B0(i)+db,D0(i)+dd)];
    % 
    %         H = [(H(:,6)+H(:,4)-2*H(:,5))./db.^2,...         % d²P/dB²
    %              (H(:,2)+H(:,8)-2*H(:,5))./dd.^2,...         % d²P/dD²
    %              (H(:,3)+H(:,7)-H(:,1)-H(:,9))./(4*db.*dd)]; % d²P/dBD
    % 
    %         % Information matrix = inv(H)
    %         H = -H./(H(:,1).*H(:,2) + H(:,3).^2); % U²D,U²B,UBD
    % 
    %         Sd(i) = 2*sqrt(H(:,1));
    %         Sb(i) = 2*sqrt(H(:,2));
    %         Sg(i) = 2*sqrt(H(:,1) + cz.^2.*(H(:,2)-2*H(:,3)));
    %     end
    % end
end