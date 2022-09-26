function [U,IAM] = uncertainty(S,x,varargin)
% U = UNCERTAINTY(S,X,[SUNAZ,SUNEL,B,T]) - Estimate expanded uncertainty [W/m²] for sensor(s) S,
%   given readings X, at solar position(s) SUNAZ,SUNEL, and a 'directional component' B, that  
%   stands for the equivalent sum of direct normal irradiance components that affect S, that is:
%
%       B = BNI + circumsolar, for GHI and GTI sensors
%       B = circumsolar, for DHI sensors
%       B = reflected BNI + reflected circumsolar, for USW sensors
%       B = BNI for BNI sensors
%
%  Expanded uncertainty is then estimated as: 
%
%  U² = ( S.zero_offset )² + 
%       ( S.calibration/100 · X )² + 
%       ( S.tilt_response/100 · X )² + 
%       ( S.temp_response/100 · | T - 25 |·X )² + 
%       ( S.cosine_response / 1000 · sin(AOI)/sin(80°) · B ) )² + 
%       ( S.azimuth_response/100 · B )²
%
% If AOI, B, or T are not provided, default values (90°, min(1000,0.9·X/cos(AOI)), 50°) will
% return a "worst-case" estimate of the effects of temperature & directional response.
%
%   TODO: currently, directional component B affects sensors also when |AOI| > 90°. While this
%       might be reasonable for sensors with a dome (glass refraction), it makes no sense for
%       flat photosensors.
%
%   TODO: generalize the sensor type 'hBNI' (pseudo-BNI, or horizontal BNI) currently detected
%       automatically only if type = BNI and model is {RSI, SPN1}.
%
% See also: METEODATA.GETUNCERTAINTY

    narginchk(2,6);
    validateattributes(x,{'numeric'},{'2d','real','size',[NaN,numel(S)]});

    ARGS = {'sunaz','sunel','B','Ts'};
    varargin(end+1:4) = {[]};
    missing = cellfun(@isempty,varargin);
    if ~all(missing)
        parsestruct(varargin(~missing),ARGS(~missing),'numeric','-r','-v','compatible',size(x));
        [varargin{~missing}] = compatiblesize(varargin{~missing});
    end
    [sunaz,sunel,B,Ts] = deal(varargin{:});
    
    % TODO: this should be a sensor type on its own!
    hBNI = contains({S.type},'BNI') & contains({S.model},{'RSI','SPN1'});
    [S(hBNI).type] = deal('hBNI');
    [S(hBNI).zero_offset] = deal(0);
    
    trashout = false;

    if isempty(sunel)
        AOI = 90;
        trashout = true;
    else
        validateattributes(sunel,{'numeric'},{'2d','real'});
        validateattributes(sunaz,{'numeric'},{'2d','real'});
        x = compatiblesize(x,sunel,sunaz);
        validateattributes(x,{'numeric'},{'size',[NaN,numel(S)]});
        AOI = sensorAOI(S,sunel,sunaz);
    end
    if isempty(Ts), Ts = 50; end
    
    x(x < 0) = 0;
    U = sqrt([S.zero_offset].^2 + 1e-4*([S.calibration].^2 + ...
        [S.tilt_response].^2 + ([S.temp_response].*abs(Ts-25)).^2).*x.^2);
    
    if size(x,1) < numel(x) || size(x,2) < numel(S)
        U = compatiblesize(U,x,1:numel(S));
    end
    
    if nargout > 1
        if nargin < 5 || isempty(B), trashout = true; end
        IAM = ones(size(U));
        assert(~trashout,'IAM factor needs solar elevation and directional fraction');
    end
    
    missingB = nargin < 5 || isempty(B);
    if missingB, B = zeros(size(U)); else, B = compatiblesize(B,U); end
    
    AOI(~isfinite(AOI)) = 90;
    assert(all(AOI >= 0,'all'))

    for j = 1:numel(S)
        if ~(S(j).azimuth_response > 0 || S(j).cosine_response > 0), return; end

        switch S(j).type
        case {'BNI','hBNI'}
            if missingB, B(:,j) = x(:,j); end
            isnormal = true;
        case {'GTI','GHI','DHI'}
            if missingB, B(:,j) = min(1000,0.9*x(:,j)./max(0,cosd(AOI(:,j)))); end
            isnormal = false;
        case {'USW'}
        % TODO: currently assuming specular reflection
            if missingB, B(:,j) = min(1000,0.9*x(:,j)./max(0,cosd(AOI(:,j)))); end
            isnormal = false;
        otherwise, error('Unknown sensor type');
        end

        if isnormal
            behind = (AOI(:,j) > 92.5);
            AOI(behind,j) = NaN; %#ok<AGROW>
            loangle = AOI(:,j) > 87.5;
            lo_aoi = AOI(loangle,j)-2.5;
            k = tand(AOI(:,j));
            k(loangle) = cosd(lo_aoi)./(1-sind(lo_aoi)); % see T20220911_low_angle_uncertainty.m
        else
            k = sind(AOI(:,j));
            % k(loangle) = 0.5/sind(87.5)*(1+sind(lo_aoi)); % difference is minimal
        end
        k = k*S(j).cosine_response/(1000*sind(80));
        k = hypot(S(j).azimuth_response/100,k);
        U(:,j) = hypot(U(:,j),B(:,j).*k);
        
        if nargout > 1 && ~isempty(S(j).fIAM)
            if isnormal
                F = 1; 
            else
                F = max(0,min(1,B(:,j).*cosd(AOI(:,j))./x(:,j))); % directional fraction
            end
            IAM(:,j) = S(j).fIAM(AOI(:,j)).*F + (1-F);
        end
    end
    if isscalar(S) && ~isscalar(U), U = reshape(U,size(x)); end
end

function AOI = sensorAOI(S,sunel,sunaz)

    assert(~any(sunel > 90),'Expecting sunel <= 90');
    if ~isscalar(S)
        AOI = zeros(numel(sunel),numel(S));
        for j = 1:numel(S)
            AOI(:,j) = sensorAOI(S(j),sunel,sunaz);
        end
        return
    end
    
    switch S.type
        case {'BNI'}, AOI = 0;
        case {'GHI','DHI','hBNI'}, AOI = 90 - sunel;
        case {'USW'}, AOI = sunel; % specular reflection
        case {'GTI'}
            c = solarposition.cos_incidence_angle(sunel,sunaz,S.mount.tilt,S.mount.azimuth);
            AOI = acosd(c);
        otherwise, error('Unknown sensor type');
    end
end
