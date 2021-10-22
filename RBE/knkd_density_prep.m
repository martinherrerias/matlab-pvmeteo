function [X0,X0_test,info] = knkd_density_prep(MD,dt)
% PROVISIONAL: generate a table with variables {cosz,kc,Kc,Kt,lastkn,lastkd} that can be used
%   to generate conditional transition PDF models for [kn,kd]. Separate this table into a test
%   and training data set.
%
% FUTURE: KNKD_DENSITY_FIT should operate directly on a MeteoData object, generating lagged/
%   averaged variables (e.g. lastkn, Kt) following some standardized syntax. e.g. 'last.kn', 
%   or function handles.
    
    AVG_WINDOW = round(1/hours(MD.timestep));
    TEST_SHARE = 0.25;

    MD = meteoQC.flagged2nan(MD,'all');
    if nargin > 1 && dt > MD.timestep
        MD.interval = 'c';
        MD = resamplestructure(MD,seconds([MD.timestep,dt]),'-centered');
    end
    if ~all(isfield(MD,{'kt','kn','kd'}))
        MD = bestimate(MD);
    end
    
    navg = @(x) circshift(movmean(x,[AVG_WINDOW-1,0],'omitnan'),1);
    % nstd = @(x) circshift(movstd(x,[AVG_WINDOW-1,0],'omitnan'),1);

    % PROVISIONAL: kn,kd now hardcoded to be the first variable pair
    X0 = struct();
    X0.kn = MD.kn;
    X0.kd = MD.kd;

    kt = MD.kn + MD.kd;
    if ~isfield(MD,'kc')
        MD = addsource(MD,'kc',kt.*MD.ENI.*sind(MD.sunel)./MD.CSGHI,strjoin(getsourceof(MD,{'kt','CSGHI'}),'/'));
        MD.kc(MD.kc < 0 | MD.kc > 1.6) = NaN;
    end

    X0.cosz = sind(MD.sunel);
    X0.Kt = navg(kt);
    % X0.sKt = nstd(kt);
    X0.kc = MD.kc;
    % X0.sKc = nstd(X0.kc);
    X0.Kc = navg(X0.kc);
    X0.lastkn = circshift(MD.kn,1);
    X0.lastkd = circshift(MD.data.kd,1);

    F = MD.flags;
    ok = cellfun(@(f) isfinite(X0.(f)),fieldnames(X0),'unif',0);
    ok = all(cat(2,ok{:}),2);
    ok = ok & all([F.kn,F.kd,F.kc,circshift(F.kn,1),circshift(F.kd,1)]== 0,2);
    ok = ok & ~MD.dark & ~MD.missing & MD.kd > 0 & (MD.flags.kt == 0);
    ok(1) = false; 
                                
    DEF = {'kn','k_n',getsourceof(MD,'kn');
           'kd','k_d',getsourceof(MD,'kd');
           'cosz','cos(\theta_z)',getsourceof(MD,'sunel');
           'Kt','\bar{k_t^H}',sprintf('movmean(kn+kd,[%d-2,-1])',AVG_WINDOW);
           'kc','k_c',getsourceof(MD,'kc');
           % 'sKc','\sigma_c^H',sprintf('movstd(kc,[%d-2,-1])',AVG_WINDOW);
           % 'sKt','\sigma_t^H',sprintf('movstd(kc,[%d-2,-1])',AVG_WINDOW);
           'lastkn','k_{n,i-1}','circshift(kn,1)';
           'lastkd','k_{d,i-1}','circshift(kd,1)'};
    DEF(cellfun(@iscell,DEF)) = [DEF{cellfun(@iscell,DEF)}];

    VAR = fieldnames(X0);
    X0 = struct2cell(X0);
    X0 = cat(2,X0{:});

    info.texlbl = parselist(VAR,DEF(:,[2,1]));
    info.descriptions = parselist(VAR,DEF(:,[3,1]));
    info.descriptions = cellfun(@(a,b) [a ' = ' b],VAR,info.descriptions,'unif',0);

    X0 = X0(ok,:);
    N = size(X0,1);

    % Divide in training/test datasets
    rng(123);
    test = false(N,1);
    test(randperm(N,ceil(N*TEST_SHARE))) = true;
    
    X0 = array2timetable(X0,'rowtimes',MD.t(ok),'variablenames',VAR);
    X0_test = X0(test,:);
    X0 = X0(~test,:);

    info.dt = MD.timestep;
    info.data = sprintf('%s (%d train / %d test points), %0.1f%% availability',...
        meteofilename(MD.location,MD.t(1),MD.t(end),MD.timestep),...
        size(X0,1),size(X0,2),nnz(ok)/nnz(~MD.dark)*100);
    
end