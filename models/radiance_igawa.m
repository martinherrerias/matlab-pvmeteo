function varargout = radiance_igawa(MD,SP,SkyRegions,model)
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
%   SunPos.Az - Nt-vector, Solar azimuth angle (degrees), convention is irrelevant, but should be
%   	coherent with the function arguments passed to LEA(az,el).
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

    if nargin == 1, model = MD;  % allow dummy call for initialization of LzEd
    elseif nargin < 4 || isempty(model), model = '2004';
    end
    assert(ischar(model) && any(strcmp(model,{'2002','2004'})),'Unknown model');
    
    % Load normalization function LzEd(Si,Gs) [inverse of integral of relative radiance] 
    % Also grid of integration points Xs with quadrature weights Ws.
    [LzEd,Xs,Ws] = getLzEd(model);
    
    if nargin <= 1, return; end
    
    narginchk(2,4);
    % tol = getSimOption('RelTol');
        
    parsestruct(MD,{'DHI','GHI','ENI'},'opt',{'AMr','kd','CSGHI'},'-n','-r','-v','-e');
    parsestruct(SP,{'El','Az'},'-n','-f','-v','size',[numel(MD.GHI),1]);
    Nt = numel(MD.GHI);
    
    % TODO: Should absolute air-mass be used??
    if ~(isfield(MD,'AMr')), MD.AMr = pvl_relativeairmass(90-SP.El); end

    useful = MD.GHI > 0 & MD.DHI > 0 & SP.El >= 0;
    dark = SP.El < 0 | MD.GHI == 0 | MD.DHI == 0;
    m = MD.AMr(useful);
    sunel = SP.El(useful);
    sunaz = SP.Az(useful);
    
%     % MHA: using pre-calculated clear-sky GHI!
%     if isfield(MD,'CSGHI')
%         Kc = min(1,MD.GHI./MD.CSGHI); 
%     else
        TL0 = 2.5;
        Ee0 = MD.ENI(useful).*max(0,cosd(sunel)); 
        Seeg = 0.84*Ee0./m.*exp(-0.027*TL0*m);      % ~ CSGHI
        Kc = MD.GHI(useful)./Seeg;                 % Clear Sky Index
        Kc = max(0,min(1,Kc));
%     end
    
    % Ce = Diffuse fraction, aka "Cloud Ratio"
    if isfield(MD,'kd'), Ce = MD.kd(useful); 
    else
        Ce = MD.DHI(useful)./MD.GHI(useful);
        Ce = max(0,min(1,Ce));
    end
   
    Ces = (m.^(0:4))*[0.01299,0.07698,-0.003857,0.0001054,-0.000001031]';
    Cle = (1-Ce)./(1-Ces);    % cloudless index
    Si = Kc+sqrt(Cle);        % Sky index
    Si = max(0,min(Si,2.0));
    
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
    clear Ces Cle Ce gs Kc m
    
    Nq = size(Xs,1);
    gg = atan2d(Xs(:,3),hypot(Xs(:,1),Xs(:,2)))'; % altitude angle of quadrature points
    
    % Classify the Nq quadrature points into the SkyRegions they fall into:
    prj = polyprojector('azim');
    prj_reg = projectstatic(SkyRegions,prj);
    prj_Xs = prj.r0*prj.fun(Xs(:,3)).*Xs(:,1:2)./hypot(Xs(:,1),Xs(:,2));

    inothers = false(Nq,1);
    B = false(Nq,SkyRegions.n.sky);
    for j = 1:numel(SkyRegions.sky)
        B(:,j) = ~inothers;
        B(~inothers,j) = insidepolygon(prj_reg.sky{j},prj_Xs(~inothers,1),prj_Xs(~inothers,2));
        inothers = inothers | B(:,j);
    end
    if SkyRegions.hasimplicit.sky && ~all(inothers)
        B(:,end) = ~inothers;
    elseif ~all(inothers)
        warning('%d points outside all regions! re-normalizing',nnz(~inothers));
        Ws = Ws.*(1+sum(Ws(~inothers))/sum(Ws));
    end
    
    D_sky = NaN(Nt,SkyRegions.n.sky);
    if ~isempty(SkyRegions.cs)
        Ncs = SkyRegions.n.solar;
        D_sun = NaN(Nt,Ncs);
    else
        D_sun = zeros(Nt,0);
    end  
    D_sky(dark,:) = 0;
    D_sun(dark,:) = 0;

    % Divide time-steps into chunks, to avoid memory issues
    Nchunks = ceil(Nt*Nq/maxarraysize('single')*4);
    ChunkSize = ceil(Nt/Nchunks);
    for k = 1:Nchunks
        idx = false(Nt,1);
        idx((k-1)*ChunkSize + 1 : min(k*ChunkSize,Nt)) = true;
        idx = useful & idx;     % index for Nt-variables
        ix = idx(useful);       % index for reduced variables

        % Angle between quadrature points and the sun, [nnz(idx),Nq] array
        zz = acosd(sph2cartV(sunaz(ix),sunel(ix))*Xs');
        
        Le = relative_radiance(Si(ix),gg,zz,model).*R(ix); % [nnz(idx),Nq] array
                
        % Calculate near-scattering (circumsolar region) components
        if ~isempty(SkyRegions.cs)
            inothers = false(size(zz));
            for j = 1:Ncs % from the smallest radius up...
                % ini = points inside ring (i), but not in ring (i-1)
                ini = zz <= SkyRegions.cs(j) & ~inothers;
                if any(ini,'all')
                    inothers = inothers | ini;
                    D_sun(idx,j) = min(Le./ini,[],2,'omitnan');
                    % D_sun(idx,j) = sum(Sc.*inj,2)./sum(inj,2);
                    % D_sun(idx,j) = sum(Sc.*sin(gg).*inj,2)./sum(sin(gg).*inj,2);
                    Le = Le - ini.*D_sun(idx,j);
                else
                    D_sun(idx,j) = 0;
                end
            end
        end

        % [~,Sc,Gr] = relative_radiance(Si,gg,zz); % [Nt,Nq] arrays
        % Sc = Sc.*R;
        % Gr = Gr.*R; % normalize
        % 
        % % Calculate near-scattering (circumsolar region) components
        % if ~isempty(SkyRegions.cs)
        %     inothers = false(size(zz));
        %     for j = 1:Ncs % from the smallest radius up...
        %         % ini = points inside ring (i), but not in ring (i-1)
        %         ini = zz <= SkyRegions.cs(j)*pi/180 & ~inothers;
        %         inothers = inothers | ini;
        %         D_sun(idx,j) = min(Sc./ini,[],2,'omitnan');
        %         % D_sun(idx,j) = sum(Sc.*inj,2)./sum(inj,2);
        %         % D_sun(idx,j) = sum(Sc.*sin(gg).*inj,2)./sum(sin(gg).*inj,2);
        %         Sc = Sc - ini.*D_sun(idx,j);
        %     end
        % else
        %     D_sun = zeros(Nt,0);
        % end
        % Le = Gr + Sc;

        for j = 1:SkyRegions.n.sky
            if ~any(B(:,j))
                D_sky(idx,j) = NaN;
            else
                D_sky(idx,j) = Le(:,B(:,j))*Ws(B(:,j))/sum(Ws(B(:,j)));
            end
            % D_sky(idx,j) = mean(Gr(idx,inj),2);
            % D_sky(idx,j) = sum(Gr(idx,inj).*sin(gg(inj)),2)./(sin(gg)*inj);
        end
    end

    % D_sky(isnan(D_sky)) = 0;
    % D_sun(isnan(D_sun)) = 0;
    varargout = {D_sky,D_sun};
end

function Le = radiance(az,el,sunel,sunaz,Si,R,model,filter)
% L = RADIANCE(AZ,EL,SUNEL,SUNAZ,Si,R,FILTER) - Return radiance values [W/m²·sr] for sky-points
%   at coordinates AZ,EL at each sun-position SUNEL,SUNAZ corresponding to sky-indices Si and 
%   normalization values R.
%
%   Meant to be passed as a function handle L = @(AZ,EL,FILTER) tied to a given set of weather
%   conditions SUNEL,SUNAZ,Si,R.

    sz = compatiblesize(sunel,sunaz,Si,R,az,el,'-size');
    assert(prod(sz) < maxarraysize('single')/4,...
        'Required %s array would exceed available memory',mat2str(sz));
    
    if nargin > 7 && ~isempty(filter)
        sunel = sunel(filter);
        sunaz = sunaz(filter);
        Si = Si(filter);
        R = R(filter);
    end

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
% [Le,Sc,Gr] = RELATIVE_RADIANCE(Si,GAMMA,ZETA) - Returns non-normalized* relative-radiance values
%   for sky-elements with elevation GAMMA and angle to the sun ZETA, for sky indices Si. 
%   The total relative radiance Le is split into a near (circumsolar) component Sc and a combined 
%   gradation + far-scattering component Gr.
%
%       Le = phi(gamma)·f(zeta) - total (non-normalized) radiance
%       Sc = phi(gamma)·(c·exp(d·zeta))
%       Gr = phi(gamma)·(1 - exp(d·pi/2) + e·cos²(zeta))
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
    
    %Sc = phi.*c.*(exp(d.*zeta) - exp(d*pi/2));          % near-scattering
    %Gr = phi.*(1 + e.*cos(zeta).^2);                    % gradation + far-scattering
    %Le = Sc+Gr; % phi(gamma)·f(zeta)
    
    varargout = {Le};
    % if nargout < 2
    %     varargout = {Le};
    % else
    %     Sc = single((1 + a.*exp(b)).*f);
    %     Gr = Le - Sc;
    %     varargout = {Le,Sc,Gr};
    % end
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