function MD = getuncertainty(MD)
% Estimate irradiance measurement uncertainties

    % Assume default sensors, if required
    MD = checksensors(MD,[],1);
    types = {MD.sensors.type};
    relevant = contains(types,MeteoData.varnames('irradiance'));

    types = unique(types(relevant));
    S = MD.sensors(relevant);

    [X,~,~,cols] = getbysource(MD,{S.ID});
    
    % TODO: tracking sensors!
    
    if all(isfield(MD,{'kn','kd','ENI'}))
        Bn = MD.data.kn.*MD.data.ENI;
        
        if isfield(MD,'F1'), CSn = MD.data.F1;  % Perez Circumsolar fraction
        else, CSn = MD.data.kn;                 % Hay-Davies approx.
        end
        CSn = CSn.*MD.data.kd.*MD.data.ENI;
        
        B = struct('BNI',Bn,'GHI',Bn+CSn,'DHI',CSn);
        if contains('GTI',types), B.GTI = B.GHI; end
        
        if contains('USW',types)
            if isfield(MD,'albedo'), B.USW = MD.data.albedo.*B.GHI; 
            else, B.USW = B.GHI;
            end
        end
    else
        B = struct('BNI',[],'GHI',[],'DHI',[],'GTI',[],'USW',[]);
    end
    
    if isempty(MD.uncertainty), MD.uncertainty = table_ish(); end
    for j = 1:numel(types)
        
        idx = find(contains({S.type},types{j}));
        U = zeros(MD.Nt,numel(idx));
        for k = idx
            U(:,cols(k)) = uncertainty(S(k),X(:,k),MD.data.sunaz,MD.data.sunel,B.(S(k).type));
        end
        MD.uncertainty.data.(types{j}) = U;
    end
end