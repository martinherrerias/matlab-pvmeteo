function MD = getsunpos(MD,check,rewrite)
% MD = getsunpos(MD) - get solar-position and -dependent variables for MeteoData MD. If some or
%   all of these variables exist, check them with METEOQC.SUNPOS, and replace any flagged values.

    if nargin < 2, check = true; end
    if nargin < 3, rewrite = false; end

    FLD = {'ENI','sunel','sunaz','hourangle','declination','TL','CSGHI','CSBNI','CSDHI','AMa'};
    dayonly = false(1,numel(FLD));
    if isfield(MD,'sunel')
       dayonly(6:10) = true;
    end

    if any(isfield(MD,FLD)) && check
    % Quality-check existing solar position data
        MD = meteoQC.sunpos(MD,'all',true);
    end

    if all(isfield(MD,FLD)) && ~rewrite
        MD = meteoQC.sunpos(MD);
        incomplete = false(MD.Nt,1);
        for j = 1:numel(FLD)
            if dayonly(j)
                incomplete = incomplete | ((MD.data.sunel > 0) & (MD.flags.data.(FLD{j}) > 0));
            else
                incomplete = incomplete | (MD.flags.data.(FLD{j}) > 0);
            end
            if size(incomplete,2) > 1, incomplete = all(incomplete,2); end
            if all(incomplete), break; end
        end
        if ~any(incomplete), MD = refresh(MD); return; end
    else
        incomplete = true(MD.Nt,1);
    end

    if isfield(MD,'Ta')
        Ta = MD.data.Ta;
        if size(Ta,2) > 1, Ta = mean(Ta,2,'omitnan'); end
        ok = meteoQC.temperature_check(Ta);
        noTa = ~ok & incomplete;
        Ta(noTa) = interp1(MD.t(ok),Ta(ok),MD.t(noTa),'nearest','extrap');
        Ta = Ta(incomplete);
        noTa = noTa(incomplete);
    else
        noTa = true(nnz(incomplete),1);
        Ta = 10;
    end

    if isfield(MD,'Patm')
        Patm = MD.data.Patm;
        if size(Patm,2) > 1, Patm = mean(Patm,2,'omitnan'); end
        ok = meteoQC.pressure_check(Patm,MD.location.altitude);
        noPa = ~ok & incomplete;
        Patm(noPa) = interp1(MD.t(ok),Patm(ok),MD.t(noPa),'nearest','extrap');
        Patm = Patm(incomplete);
        noPa = noPa(incomplete);
    else
        Patm = pvl_alt2pres(MD.location.altitude);
        noPa = true(nnz(incomplete),1);
    end

    if ~isfield(MD,'ENI')
        MD = addsource(MD,'ENI',solarposition.extraradiation(MD.t),'sp.extraradiation');
    end

    if isfield(MD,'TL')
        TL = MD.data.TL;
        bad = MD.flags.data.TL > 0 & incomplete;
        if any(bad)
            TL(bad) = linketurbidity(MD.location,MD.t(bad));
            % MD.flags = checkfield(MD.flags,MD,@(x) ~bad,'NA','TL',true);
        end
        TL = TL(incomplete);
    else
        TL = linketurbidity(MD.location,MD.t(incomplete));
    end

    [SP,C] = solarposition.effectivesunpos(MD.location,MD.t(incomplete),Patm,Ta,TL,...
        'step',MD.timestep,'interval',MD.interval,'N2E');           
    SP = renamefields(SP,{'Az','sunaz';'El','sunel';'w','hourangle';'dec','declination'});
    SP = completestruct(SP,C);
    SP.TL = TL;

    SP.AMa = NaN(size(SP.sunel));
    SP.AMa(SP.sunel >= 0) = pvl_relativeairmass(90-SP.sunel(SP.sunel >= 0));
    SP.AMa = SP.AMa.*Patm/101325;

    if (MD.interval == 'i' || MD.timestep < minutes(5)) && any((noTa | noPa) & ~(SP.sunel < -1))
    % Estimate the error due to refraction correction using Patm/Ta defaults:
    %  + up to ~14% error in barometric pressure
    %  + up to ~12% error due to temperature (-20°C instead of 10°C)
    % ... and flag those points (sunel/sunaz) for which the error

        MIN_DW = 1e-3; % ~1/4 second (SPA precision is ~3e-4)

        if MD.interval == 'i', dw = 0;
        else, dw = 360*days(MD.timestep)/4;
        end
        dw = max(dw,MIN_DW);

        r = solarposition.refraction(SP.sunel,'-app');
        r = (0.14.*noPa + 0.12.*noTa).*r; 
        toflag = (r > dw) & ~(SP.sunel < -1);
    else
        toflag = false;
    end

    [b,MD.flags] = flagbit(MD.flags,'interp');
    for fld = fieldnames(SP)'
        f = fld{1};
        switch f
        case {'sunel','sunaz','hourangle','declination','ENI'}, tag = ['effective.' f];
        case {'CSGHI','CSBNI','CSDHI'}, tag = ['ineichen.' f];
        case 'AMa', tag = 'kastenyoung1989';
        case 'TL', tag = 'linketurbidity';
        end

        if ~isfield(MD,f) || rewrite
            MD = addsource(MD,f,NaN(MD.Nt,1),tag,true);
            MD.data.(f) = SP.(f);
            tochange = true(MD.Nt,1);
            ic = 1;
        else
        % Don't modify existing columns, except if they use the same model

            c = getsourceof(MD,f);
            ic = contains(c,tag);
            
            
            tochange = false(MD.Nt,1);
            if any(ic)
                tochange(incomplete) = MD.flags.data.(f)(incomplete,ic) ~= 0 | ...
                      ~isfinite(MD.data.(f)(incomplete,ic));
                MD.data.(f)(tochange,ic) = SP.(f)(tochange(incomplete));
            end
        end
        if contains(f,{'sunel','sunaz','AMa'})
            MD.flags.data.(f)(tochange,ic) = bitset(MD.flags.data.(f)(tochange,ic),b,toflag);
        end

        if all(tochange)
            MD = setsourceof(MD,f,tag);
        else
            c = getsourceof(MD,f);
            if ~contains(c{1},tag)
                MD = setsourceof(MD,f,strjoin([c,{tag}],'/'));
            end
        end
    end
    MD = refresh(MD);
end