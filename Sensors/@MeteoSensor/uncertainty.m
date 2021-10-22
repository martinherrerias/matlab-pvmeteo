function U = uncertainty(S,x,varargin)
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

    if isempty(sunel), AOI = 90;
    else
        x = compatiblesize(x,sunel,sunaz);
        AOI = sensorAOI(S,sunel,sunaz);
    end
    if isempty(Ts), Ts = 50; end
    
    x(x < 0) = 0;
    U = sqrt([S.zero_offset].^2 + ([S.calibration]/100.^2 + ...
        [S.tilt_response]/100.^2 + ([S.temp_response]/100.*abs(Ts-25)).^2).*x.^2);
    
    if size(x,1) < numel(x) || size(x,2) < numel(S)
        U = compatiblesize(U,x,1:numel(S));
    end

    for j = 1:numel(S)
        if ~(S(j).azimuth_response > 0 || S(j).cosine_response > 0), return; end

        if nargin < 5 || isempty(B)
            switch S(j).type
                case {'BNI'}, B = x;
                case {'GTI','GHI','DHI','USW'}, B = min(1000,0.9*x./max(0,cosd(AOI(:,j)))); 
            otherwise, error('Unknown sensor type');
            end
        end
        B = compatiblesize(B,U);
        
        k = min(1,hypot(S(j).azimuth_response/100,S(j).cosine_response/1000.*sind(AOI(:,j))/sind(80)));
        U(:,j) = hypot(U(:,j),B(:,j).*k);
    end
    if isscalar(S) && ~isscalar(U), U = reshape(U,size(x)); end
end

function AOI = sensorAOI(S,sunel,sunaz)

    if ~isscalar(S)
        AOI = zeros(numel(sunel),numel(S));
        for j = 1:numel(S)
            AOI(:,j) = sensorAOI(S(j),sunel,sunaz);
        end
        return
    end
    
    switch S.type
        case {'BNI'}, AOI = 0;
        case {'GHI','DHI'}, AOI = 90 - sunel;
        case {'USW'}, AOI = - sunel;
        case {'GTI'}
            c = solarposition.cos_incidence_angle(sunel,sunaz,S.mount.tilt,S.mount.azimuth);
            AOI = acosd(c);
        otherwise, error('Unknown sensor type');
    end
end
