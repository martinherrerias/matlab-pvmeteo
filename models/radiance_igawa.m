function varargout = radiance_igawa(MD,SP,SkyRegions,model,ikb)
% LEA = RADIANCE_IGAWA(METEODATA,SUNPOS) - Estimate diffuse-radiance distribution function
%    LEA(az,el) for the sky conditions and time-steps defined by METEODATA and SUNPOS, using
%    the Igawa et al. (2004) radiance distribution model.
%
% [DSKY,DSUN] = RADIANCE_IGAWA(METEODATA,SUNPOS,SKYREGIONS) - Estimate irradiance-component weight
%    factors (i.e. mean sky-region radiance values), for an arbitrary sky-tessellation SKYREGIONS
%    at each of the time-steps defined by METEODATA and SUNPOS, using the Igawa et al. (2004) all-
%    weather radiance distribution model.
%
%    With this syntax, RADIANCE_IGAWA uses SPHQUADRATUREPOINTS to evaluate N ~ 1/SimOptions.RelTol
%    quasi-uniformly distributed points on the hemispherical unit dome. Points are classified into
%    the different SKYREGIONS, and the results DSKY correspond to the weighted radiances of points
%    in each given region.
%
%    The extraction of a circumsolar component DSUN is somewhat artificial. Radiance in the model
%    is the _product_ of a gradation function (on sky-element elevation) and a scattering function
%    (on sky-element angle vs the sun), and cannot be cleanly split into a _sum_ of components.
%    DSUN is therefore set as the MINIMUM radiance within any given circumsolar ring, so that it
%    can be substracted from the remaining DSKY (background) radiance without risk of producing
%    non-physical (i.e. negative) values.
%
%    NOTE: the use of (a) circumsolar-region(s) is strongly discouraged! It requires a much higher
%    computational burden, as it must be calculated on a time-step basis, and for most sky-
%    tessellations it actually increases projection- and averaging errors. Only for very coarse 
%    sky-tessellations (e.g. 'isotropic', 'sunto',...), a 25-30° radius does (slightly) improve
%    projection accuracy.
%
% INPUT:
% 	MeteoData.GHI - Nt-vector, Global Horizontal Irradiance [kWh/m²]
%	MeteoData.DHI - Nt-vector, Diffuse Horizontal Irradiance [kWh/m²]
%	MeteoData.ENI - Nt-vector, Normal Extraterrestrial irradiance [kWh/m²]
%	MeteoData.AMr - (optional) Nt-vector, Relative Air Mass
%
%	SunPos.El - Nt-vector, Solar elevation angle (degrees)
%   SunPos.Az - Nt-vector, Solar azimuth angle (degrees), convention is assumed to be N2E 
%       (astronomical), for SKYREGIONS defined in a system with x = East, y = North.
%       For direct calls to LEA(az,el), convention is irrelevant, but should be consistent.
%
% OUTPUT:
%   LEA = @(AZ,EL,[FILTER]) takes 1·M horizontal vectors of sky-element azimuth and elevation angles  
%       AZ, EL [degrees] and returs an array of Nt·M absolute diffuse-radiance values [kWh/m²·sr], 
%       for the M sky elements in AZ,EL and Nt time-steps in METEODATA/SUNPOS.
%
%       NOTE that the provided function is bound to the provided METEODATA/SUNPOS. The optional 
%       third argument FILTER can be used to calculate & return only a subset of rows (i.e. only
%       the timesteps in FILTER). This is useful when a large number of sky-elements AZ,EL need 
%       to be calculted simultaneously, but doing so for all Nt time-steps would require an array
%       that exceeds the available memory. For repeated calls using the same FILTER, consider 
%       using instead:
%
%           MD = FILTERSTRUCTURE(METEODATA,FILTER);
%           SP = FILTERSTRUCTURE(SUNPOS,FILTER);
%           LEA = RADIANCE_IGAWA(MD,SP); 
%
%   AUX.f - @(zeta) scattering-indicatrix function, for angle from sun zeta (rad)
%   AUX.phi - @(gamma) gradation function for elevation angle gamma (rad)
%   AUX.LzEd - normalization values (inverse of integration over sky-dome).
%   AUX.Kc, .Ces, .Cle, .Si - Clear sky index, std. cloud ratio, cloudless index, sky index.
%
%   DSKY,DSUN - [Nt,SKYREGIONS.n.sky] and [Nt,SKYREGIONS.n.solar] arrays of mean diffuse-radiance
%       components [kWh/m²·sr] for polygonal regions SKYREGIONS.sky and SKYREGIONS.solar.
%
% N. Igawa, Y. Koga, T. Matsuzawa, and H. Nakamura, “Models of sky radiance distribution and sky
% luminance distribution,” Solar Energy, vol. 77, no. 2, pp. 137–157, 2004.
%
% EXAMPLE:
%     [MD,T,Loc] = getMeteoData();
%     [MD,T,SP,Loc] = completemeteodata(MD,T,Loc,'filter','12days');
%     SR = ShadingRegions('rb32');
%     [Dsky,Dsun] = radiance_igawa(MD,SP,SR);
%
% See also: PVLMOD_HAYDAVIES, POAIRRADIANCE

    if nargin == 0, test(); return; end

    if nargin == 1 && ischar(MD), model = MD;  % allow dummy call for initialization of LzEd
    elseif nargin < 4 || isempty(model), model = '2004';
    end
    if nargin < 5, ikb = false; end
    assert(ischar(model) && any(strcmp(model,{'2002','2004'})),'Unknown model');
    
    % Load normalization function LzEd(Si,Gs) [inverse of integral of relative radiance] 
    % Also grid of integration points Xs with quadrature weights Ws.
    [LzEd,Xs,Ws] = getLzEd(model);
    Nq = numel(Ws);
    
    if isa(MD,'MeteoData')
        [MD,SP] = MD.legacy();
        narginchk(1,5);
    else
        if nargin <= 1, return; end
        narginchk(2,5);
    end
    % tol = getSimOption('RelTol');
 
    parsestruct(MD,{'DHI','GHI','ENI'},'opt',{'AMr','CSGHI','CSDHI'},'-n','-r','-v','-e');
    parsestruct(SP,{'El','Az'},'-n','-f','-v','size',[numel(MD.GHI),1]);
    Nt = numel(MD.GHI);

    useful = MD.GHI > 0 & MD.DHI > 0 & SP.El >= 0;
    dark = SP.El < 0 | MD.GHI <= 0 | MD.DHI <= 0;
    sunel = single(SP.El(useful));
    sunaz = single(90 - SP.Az(useful)); % E2N
    
    Si = skyindex(MD,SP,useful,ikb);
    
    R = MD.DHI(useful).*LzEd(Si,sunel)./relative_radiance(Si,90,90-sunel,model);
    
    if nargin < 3 || isempty(SkyRegions)
        R = revertfilter(R,useful,[],NaN); R(dark) = 0;
        Si = revertfilter(Si,useful,[],NaN); Si(dark) = 0;
        LEA = @(az,el,varargin) radiance(az,el,SP.El,SP.Az,Si,R,model,varargin{:});
        varargout = {LEA};
        if nargout > 1
            AUX = struct('Si',Si,'LzEd',LzEd(Si,SP.El));
            [AUX.phi,AUX.f] = relative_radiance(Si,model);
            varargout{2} = AUX;
        end
        return
    end
    
    D_sky = NaN(Nt,SkyRegions.n.sky);
    if ~isempty(SkyRegions.cs)
        D_sun = NaN(Nt,SkyRegions.n.solar);
    else
        D_sun = zeros(Nt,0);
    end  
    D_sky(dark,:) = 0;
    D_sun(dark,:) = 0;
    
    % Initialize DISCRETE_RADIANCE, return altitude angle of quadrature points
    gg = discrete_radiance(Xs,Ws,SkyRegions);

    % Divide time-steps into chunks, to avoid memory issues
    Nt = nnz(useful);
    Nchunks = ceil(Nt*Nq/maxarraysize('single')*8);
    ChunkSize = ceil(Nt/Nchunks);
    for k = 1:Nchunks
        
        ix = false(Nt,1);
        ix((k-1)*ChunkSize + 1 : min(k*ChunkSize,Nt)) = true;
        idx = useful;
        idx(useful) = ix; % index for Nt-variables
        
%         idx = false(Nt,1);
%         idx((k-1)*ChunkSize + 1 : min(k*ChunkSize,Nt)) = true;
%         idx = useful & idx;     % index for Nt-variables
%         ix = idx(useful);       % index for reduced variables

        % Angle between quadrature points and the sun, [nnz(idx),Nq] array
        zz = acosd(sph2cartV(sunaz(ix),sunel(ix))*Xs');
        
        % Radiance at quadrature points
    	Le = relative_radiance(Si(ix),gg',zz,model).*R(ix); % [nnz(idx),Nq] array

        % Discretized radiance for each element of SkyRegions
        [D_sky(idx,:),D_sun(idx,:)] = discrete_radiance(Le,zz);
    end

    % D_sky(isnan(D_sky)) = 0;
    % D_sun(isnan(D_sun)) = 0;
    varargout = {D_sky,D_sun};
end

function Si = skyindex(MD,SP,useful,ikb)

    % TODO: Should absolute air-mass be used??
    if ~(isfield(MD,'AMr')), MD.AMr = pvl_relativeairmass(90-SP.El); end

    m = MD.AMr(useful);
    
    % Clear Sky Index
    if ikb && isfield(MD,'CSGHI')
        Kc = MD.GHI(useful)./MD.CSGHI(useful); 
    else
        Kc = MD.GHI(useful)./approxCSGHI(MD.ENI(useful),m); 
    end
    Kc = max(0,min(1,Kc));
    
    % Ce = Diffuse fraction, aka "Cloud Ratio"
    Ce = MD.DHI(useful)./MD.GHI(useful);
    Ce = min(Ce,1);
   
    if ikb && all(isfield(MD,{'CSDHI','CSGHI'}))
        Ces = MD.CSDHI(useful)./MD.CSGHI(useful); 
    else    
        % MHA: Igawa et al. use only points with sunel > 5°, so m < ~10.3 and Ces < 0.5
        m = min(20,m);
        
        Ces = (m.^(0:4))*[0.01299,0.07698,-0.003857,0.0001054,-0.000001031]';
    end
    Cle = (1-Ce)./(1-Ces);    % cloudless index
    Si = Kc+sqrt(Cle);        % Sky index
    Si = max(0,min(single(Si),2.0));
end

function Seeg = approxCSGHI(Ee0,m)
    TL0 = 2.5;
    Seeg = 0.84*Ee0./m.*exp(-0.027*TL0*m);  % ~ CSGHI
end

function Le = radiance(az,el,sunel,sunaz,Si,R,model,filter)
% L = RADIANCE(AZ,EL,SUNEL,SUNAZ,Si,R,FILTER) - Return radiance values [W/m²·sr] for sky-points
%   at coordinates AZ,EL at each sun-position SUNEL,SUNAZ corresponding to sky-indices Si and 
%   normalization values R.
%
%   Meant to be passed as a function handle L = @(AZ,EL,FILTER) tied to a given set of weather
%   conditions SUNEL,SUNAZ,Si,R.

    if nargin > 7 && ~isempty(filter)
        sunel = sunel(filter);
        sunaz = sunaz(filter);
        Si = Si(filter);
        R = R(filter);
    end
    
    sz = compatiblesize(sunel,sunaz,Si,R,az,el,'-size');
    assert(prod(sz) < maxarraysize('single')/4,...
        'Required %s array would exceed available memory',mat2str(sz));

    % Get angle from elements to sun
    zeta = sph2cartV(sunaz,sunel)*sph2cartV(az,el)';
    zeta = acosd(zeta);

    % [Le,Sc,Gr] = relative_radiance(Si,el,zeta); 
    % Le = Le.*R;
    % Sc = Sc.*R;
    % Gr = Gr.*R;
    
    Le = relative_radiance(Si,el,zeta,model).*R;
    Le(:,el < 0) = 0;
end

function varargout = relative_radiance(Si,varargin)
% [PHI,F] = RELATIVE_RADIANCE(Si) - Returns relative gradation and scattering-indicatrix functions
%   for sky-index/indices Si. The gradation function PHI = @(GAMMA) where GAMMA is the altitude of
%   a sky element (degrees), expresses changes of luminance along the sky vertical meridian caused
%   by diffusion. The scattering indicatrix F = @(ZETA) where ZETA is the angle between a sky
%   element and the sun (degrees), represents the distribution of scattered sun beams diffused 
%   into the surrounding space of an atmospheric particle. The radiance of a sky-element (g,z) at
%   the given Si is PROPORTIONAL* to PHI(g)·F(z).
%
% Le = RELATIVE_RADIANCE(Si,GAMMA,ZETA) - Returns non-normalized* relative-radiance values
%   for sky-elements with elevation GAMMA and angle to the sun ZETA, for sky indices Si. 
%
%       Le = phi(gamma)·f(zeta) - total (non-normalized) radiance
%
%   (*) NOTE: values are not normalized! to get absolute radiance, they must be normalized to the
%   calculated relative radiance at zenith RELATIVE_RADIANCE(Si,90,90 - Gs) for solar elevation Gs
%   and then multiplied by the actual radiance at zenith LeZ(Si,Gs) = DHI·LzEd(Si,Gs). That is:
%
%       L = DHI·LzEd(Si,Gs)·RELATIVE_RADIANCE(Si,GAMMA,ZETA)/RELATIVE_RADIANCE(Si,90,90 - Gs)
%
%   The normalization factor LzEd(Si,Gs) is calculated by numerically integrating PHI(g)·F(z)
%   over the complete sky-dome (see GETLZED).
%
% N. Igawa, Y. Koga, T. Matsuzawa, and H. Nakamura, “Models of sky radiance distribution and sky
% luminance distribution,” Solar Energy, vol. 77, no. 2, pp. 137–157, 2004.

    narginchk(2,4);
    model = varargin{end};
    assert(ischar(model) && any(strcmp(model,{'2002','2004'})),'Unknown model');
    
    switch model
    case '2002'
        
        Si = max(0,min(2.2,Si)); % for real c, allow Si < 2.3

        % Coefficients in Igawa et al. 2002
        a = 9.25./(1+0.29*exp(3.51*Si))-1.16;
        b = 0.927./(1+8.94*exp(-2*Si))-1.16;
        c = 1.7*(0.94*Si).^3.2.*exp(1.04*Si).*(2.3-Si).^1.63;
        d = -3.5./(1+8.6*exp(-2.3*Si));
        e = 0.56./(1+89*exp(-3.95*Si));
            
    case '2004'

        Si = max(0,min(2.0,Si)); % for real c, allow Si < 2.1
    
        % Coefficients in Igawa et al. 2004
        a = 4.5./(1+0.15*exp(3.4*Si))-1.04;
        b = -1./(1+0.17*exp(1.3*Si))-0.05;
        c = 1.77*(1.22*Si).^3.56.*exp(0.2*Si).*(2.1-Si).^0.8;
        d = -3.05./(1+10.6*exp(-3.4*Si));
        e = 0.48./(1+245*exp(-4.13*Si));
    end

    if nargin < 4
        phi = @(g) (1 + a.*exp(b./sind(max(0,g))))./(1 + a.*exp(b));   % relative gradation
        f = @(z) 1 + c.*(exp(d.*z*pi/180)-exp(d*pi/2))+e.*cosd(z).^2;  % scattering indicatrix
        varargout = {phi,f};
        return;
    else
       [gamma,zeta] = deal(varargin{1:2}); 
    end
    
    phi = (1 + a.*exp(b./sind(max(0,gamma))))./(1 + a.*exp(b));     % relative gradation
    f = 1 + c.*(exp(d.*zeta*pi/180)-exp(d*pi/2))+e.*cosd(zeta).^2;  % scattering indicatrix
    Le = single(phi.*f);
    varargout = {Le};
end

function [LzEd,Xs,Ws] = getLzEd(model)
% Return an interpolant function LzEd(Si,Gs) = Lez(Si,Gs)/Ed = 1/Ih(Si,Gs) to normalize the result
% of relative_radiance Le(Si,g,z). Lez = zenith radiance, Ed = horizontal diffuse irradiance, and 
% Ih is the double integral of the relative radiance: Ih(Si,Gs) = §[Le(Si,g,z)·sin(el)·dA]
% over the complete sky dome, for a grid of solar elevation angles Gs and Sky-index Si. 

    TOL = getSimOption('RelTol');

    file = sprintf('igawa_%s_LzEd_%0.0e.mat~',model,TOL);
    file = fullfile(fileparts(mfilename('fullpath')),file);
    if ~isempty(dir(file))
        load(file,'-mat','LzEd','Xs','Ws');
        return;
    end

    fprintf('Initializing radiance_igawa...\n');
        
    [Xs,Ws] = sphquadraturepoints(1/TOL);
    assert(abs(sum(Ws)/(2*pi)-1) < TOL);
    
    gg = atan2d(Xs(:,3),hypot(Xs(:,1),Xs(:,2)));
    Nq = numel(gg);
    
    switch model
        case '2002', Si0 = single(0:0.2:2.2);
        case '2004', Si0 = single(0:0.2:2.0);
    end
    Gs0 = single(0:15:90);  
    
    [x,V] = fitinterpolant(@F,{Si0,Gs0},TOL,'makima');
    LzEd = griddedInterpolant(x,V,'makima');

    function LzEd = F(Si,Gs)
    % Integrate relative_radiance(Si,..) over the sky-dome, using quadrature points Xs
    
        n = numel(Gs);
        Ih = zeros([n,1],'single');
        Nchunks = ceil(Nq*n/maxarraysize('single')*4);
        ChunkSize = ceil(Nq/Nchunks);
        for j = 1:Nchunks
            idx = false(Nq,1);
            idx((j-1)*ChunkSize + 1 : min(j*ChunkSize,Nq)) = true;
            zz = acosd(Xs(idx,:)*[cosd(Gs),zeros(size(Gs)),sind(Gs)]');
            Le = relative_radiance(Si',gg(idx),zz,model)./relative_radiance(Si',90,90-Gs',model);
            Ih = Ih + Le'*(Xs(idx,3).*Ws(idx));
        end
        LzEd = 1./Ih;
    end
    save(file,'LzEd','Xs','Ws');
    
    % contourf(LzEd.GridVectors{1},LzEd.GridVectors{2}*180/pi,LzEd.Values'); colorbar();
end

function varargout = discrete_radiance(varargin)
% G = DISCRETE_RADIANCE(XS,WS,SR) - initialization call.
% [DSKY,DSUN] = DISCRETE_RADIANCE(LE,Z) - Estimate integrated radiance LE over the discrete regions
%   SR.sky, and SR.solar of SR. Requires a previous initialization call (above).
%
% INPUT:
%   XS,WS - quadrature points and weights, as returned by GETLZED
%   SR - SHADINGREGIONS object
%   LE - [N,M] array of radiance values, already normalized!, M must match size(XS,1)
%   Z - [N,M] array of angles between quadrature points Xs and the sun
%     G = DISCRETE_RADIANCE(Xs,Ws,SR);
%     Z = acosd(sph2cartV(sunaz,sunel)*Xs');
%     LE = relative_radiance(Si,G,Z,model).*R;
%     [DSKY,DSUN] = DISCRETE_RADIANCE(LE,Z);
        
    persistent Xs
    persistent Ws
    persistent SR
    persistent B
    
    if nargin == 3
        [Xs,Ws,SR] = deal(varargin{:});
        validateattributes(Ws,'numeric',{'real','size',[NaN,1]});
        validateattributes(Xs,'numeric',{'real','size',[numel(Ws),3]});

        [B,Ws] = assignpoints(SR,Xs,Ws);
        varargout{1} = atan2d(Xs(:,3),hypot(Xs(:,1),Xs(:,2)));
        return;
    else
        [Le,zz] = deal(varargin{:});
        Nt = size(Le,1);
        validateattributes(Le,'numeric',{'real','size',[NaN,numel(Ws)]});
        validateattributes(zz,'numeric',{'real','size',[Nt,numel(Ws)]});
    end
    
    D_sky = NaN(Nt,SR.n.sky);
    if ~isempty(SR.cs)
        D_sun = NaN(Nt,SR.n.solar);
    else
        D_sun = zeros(Nt,0);
    end  
  
    if isempty(SR.cs)
    % Calculate sky components as weighted averages of quadrature points within each region.
    
        for j = 1:SR.n.sky
            inj = B(:,j);
            D_sky(:,j) = Le(:,inj)*Ws(inj)/sum(Ws(inj));
        end
    else
    % For overlapping circumsolar regions, the problem becomes: 
    %
    %   [D_sky(t,:), D_sun(t,:)]' = x(t) = argmin( [B C(t)] W x = Le(t,:)' )
    %
    % Where B(i,j) is a boolean matrix assigning quadrature points Xs to sky regions,
    % C(t) would be a matrix (similar to B) for circumsolar regions, at time t,
    % W = diag(Ws)/mean(Ws) applies weights to each quadrature point,
    % and Le(t,:)' is the vector of radiances calculated at Xs for time t.
    %
    % The least squares solution x = [B'W²B B'W²C; C'W²B C'W²C]\[B C]'W Le is not guaranteed to
    % be physical (radiance components might be negative). To avoid this issue, let instead:
    %
    %   x(t) = [D_sky - min(Le), D_sun]', x(t) = argmin( [B C] W x = Le - min(Le) ) | x > 0
    %
    % NOTE: there is no restriction on circumsolar irradiance over several rings being monotonic!

        Ws = Ws/mean(Ws);
        
        Le_min = min(Le,[],2);
        Le = Le - Le_min;
        
        BC = zeros(Nt,SR.n.sky,SR.n.solar,'double');
        
        Y = zeros(Nt,SR.n.sky+SR.n.solar,'double');
        Y(:,1:SR.n.sky) = (Le.*Ws')*B;
        
        CC = zeros(Nt,SR.n.solar,'double');
        inothers = false(size(zz));
        for j = 1:SR.n.solar
            inj = (zz <= SR.cs(j) & ~inothers);
            inothers(inj) = true;
            wj = inj.*Ws';
            CC(:,j) = sum(wj.^2,2);
            Y(:,j+SR.n.sky) = sum(wj.*Le,2);
            for k = 1:SR.n.sky
                BC(:,k,j) = sum(wj(:,B(:,k)).^2,2);
            end
        end

        BC = permute(BC,[2,3,1]);
        BB = double(diag(sum(B.*Ws.^2,1))); % = (B.*Ws)'*(B.*Ws);

        X = arrayfun(@(j) lsqnonneg([BB,BC(:,:,j);BC(:,:,j)',diag(CC(j,:))],Y(j,:)'),1:Nt,'unif',0);
        X = cat(2,X{:});
        
        D_sky = X(1:SR.n.sky,:)' + Le_min;
        D_sun = X(SR.n.sky+1:end,:)';
    end

    varargout = {D_sky,D_sun};
end

function [B,Ws] = assignpoints(SR,Xs,Ws)
% Classify quadrature points Xs into the Sky-Regions they fall into, return a boolean matrix B,
% where B(i,j) = true implies point Xs(i,:) is inside SR.sky{j}.

    Nq = size(Xs,1);

    prj = polyprojector('azim');
    prj_reg = projectstatic(SR,prj);
    prj_Xs = prj.r0*prj.fun(Xs(:,3)).*Xs(:,1:2)./hypot(Xs(:,1),Xs(:,2));

    inothers = (Ws <= 0);
    B = false(Nq,SR.n.sky);
    for j = 1:numel(SR.sky)
        B(:,j) = ~inothers;
        B(~inothers,j) = insidepolygon(prj_reg.sky{j},prj_Xs(~inothers,1),prj_Xs(~inothers,2));
        inothers = inothers | B(:,j);
    end
    if SR.hasimplicit.sky && ~all(inothers)
        B(:,end) = ~inothers;
    elseif ~all(inothers)
        warning('%d points outside all regions! re-normalizing',nnz(~inothers));
        Ws = Ws.*(1+sum(Ws(~inothers))/sum(Ws));
    end
end

function test()

    SR = ShadingRegions('sat13_CS_10_20_30');
    LzEd = getLzEd('2004');
    N = 512;
        
    prj = polyprojector('ortographic','tol',1e-6);
    [xx,yy] = meshgrid(linspace(-1,1,N),linspace(-1,1,N));
    V = prj.inv([xx(:),yy(:)]')';
    outside = ~(V(:,3) > max(getSimOption('RelTol'),prj.cX)); 
    V = V(~outside,:);
    
    gg = discrete_radiance(V,ones(nnz(~outside),1)/nnz(~outside),SR);

    fh = GUIfigure('radiance_igawa_test','-silent');
    try close(fh); end
    
    GUIfigure('radiance_igawa_test','radiance_igawa test','3:1'); clf();
    ax = arrayfun(@(j) subplot(1,3,j),1:3);
    
    himg = imagesc(ax(1),xx(1,1:N),yy(1:N,1),ones(N),'AlphaData',reshape(~outside,[N,N]));
    axis(ax(1),'equal');
    set(ax(1),'ydir','normal');
    axis(ax(1),[-1,1,-1,1]);
    % set(ax(1),'xdir','reverse');
    colorbar(ax(1));
    
    ax(1).Visible = 'off';
    ax(2).Visible = 'off';
    ax(2).Position = ax(2).Position*[1 0 0 0; 0 1 0 0.1; 0 0 1 0; 0 0 0 0.8]';
    
    opt.rotation = eye(3);
    opt.prj = prj;
    opt.ax = ax(3);
    [~,hsky] = plot(SR,opt);
    % set(ax(3),'xdir','reverse');
    ax(3).Position(3:4) = ax(1).Position(3:4);
    
    colorlist = colormap(ax(1));
    ncolors = size(colorlist,1);
    cfun = @(x,b) colorlist(round((x-b(1))/diff(b)*(ncolors-1)+1),:);
    
    CTL = plotcontrols({'popupmenu','slider','slider','slider','text'},...
                 {'model','sunel','sunaz','Si','output'},...
                 {{'2004','2002'},[0,90],[-180,180],[0,2.1],[]},...
                 {'2004',45,120,1.2,'foo'},@update,...
                 'position',ax(2).Position,'-skipupdate','-smooth');
             
    LzEd = [];
    Le = [];
    augmented = {};
    hcs = {};
    
    update('2004',45,120,1.2)

    function update(model,sunel,sunaz,Si,~,caller,~)
        if isempty(LzEd) || strcmp(caller.Style,'popupmenu')
            [LzEd,~] = getLzEd(model);
            switch model
                case '2002', CTL(4).Max = 2.2;
                case '2004', CTL(4).Max = 2.0; CTL(4).Value = min(CTL(4).Value,2.0); 
            end
        end
        
        zz = acosd(V*sph2cartV(90-sunaz,sunel)');
        Le = LzEd(Si,sunel)*relative_radiance(Si,gg,zz,model)./relative_radiance(Si,90,90-sunel',model);
        himg.CData(~outside) = Le;
        
        [Dsky,Dcs] = discrete_radiance(Le',zz');
        
        CTL(5).String = sprintf('sum = %0.3f',nanmean(Le)*pi);
        
        clim = [min(min(Le),min(Dsky)),max(max(Le),max(Dsky))];
        
        caxis(ax(1),clim);
        opt.colors.sky = cfun(Dsky',clim);
        opt.sunpos = [sunaz,sunel];
            
        needsredraw = nargin < 6 || contains(caller.UserData.label,{'sunel','sunaz'});
        if needsredraw
            cla(ax(3));
            [~,hsky] = plot(SR,opt);
            
            % Intersect circumsolar and sky regions to color according to augmented radiance
            if SR.n.solar > 0
                sky = polygon3d.vf2poly(hsky.sky.Vertices,hsky.sky.Faces);
                cs = polygon3d.vf2poly(hsky.solar.Vertices,hsky.solar.Faces);
                ncs = numel(cs);
                augmented = cell(ncs,1);
                hcs = cell(ncs,1);
                for j = numel(ncs):-1:1
                    p = arrayfun(@(p) intersectpolygons(cs(j),p),sky,'unif',0);
                    augmented{j} = ~cellfun(@isempty,p);
                    [v,~,p] = poly2vef(cat(1,p{augmented{j}}),1);
                    hcs{j} = patch('faces',p,'vertices',v,'EdgeColor','w',...
                            'FaceColor','flat','FaceAlpha',1,'LineWidth',0.1);
                end
            end    
        end
        
        if any(Dcs + Dsky' > clim(2),'all')
            clim(2) = max(Dcs + Dsky',[],'all');
            caxis(ax(1),clim);
            opt.colors.sky = cfun(Dsky',clim);
        end
        hsky.sky.FaceVertexCData = opt.colors.sky;
        if SR.n.solar > 0
            for j = 1:numel(hcs)
                cscolors = cfun(Dcs(j) + Dsky(augmented{j}),clim);
                hcs{j}.FaceVertexCData = cscolors;
            end
        end    
    end
end