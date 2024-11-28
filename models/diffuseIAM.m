function varargout = diffuseIAM(fIAM,tilt,varargin)
% [ISO,ALB,HB,CS] = DIFFUSEIAM(FIAM) - Integrate IAM function under an isotropic radiance & flat
%   terrain assumption, returning function handles for the tilted components:
%
% INPUT:
%   FIAM - function handle, parsed by CHECKIAM.
% OUTPUT:
%   ISO(T) - Isotropic component IAM as a function of surface tilt T. Use ISO(T+A) to consider
%       infinite-row shading with shade-free angle A.
%   ALB(T) - Albedo component IAM as a function of surface tilt T
%   HB(T) - Horizon-brightening component IAM as a function of surface tilt T
%   CS(IA,R) - Approx. circumsolar component IAM, as a function of incidence angle IA, and
%       and circumsolar radius R.*
%
%       (*) NOTE: shading by the horizon is not considered for CS, only the effect of clipping
%           by the surface plane.
%
% [ISO,ALB,HB,CS] = DIFFUSEIAM(FIAM,TILT,IA,R) - return values for specific surface TILT, solar
%   incidence angle IA, and circumsolar radius R.
%
% INPUT: 
%   TILT,IA - singleton-expansion-compatible arrays of surface tilt angle(s) and solar incident 
%       angle(s). The latter can be omitted if CS is not required. 
%   R - scalar, circumsolar radius. Default = 25.
%
% OUTPUT: ISO, ALB, HB, CS will be numerical arrays, evaluated at the specific angles.
%
% TODO: integral2 crashes with discontinuous functions, e.g. checkIAM('hcpv')
%
% See also: CHECKIAM, CIRCUMSOLARIAM

    DZ = 1.0; % degrees, resolution of interpolation structure for f.iso,...
    
    fIAM = checkIAM(fIAM);

    df = @(x,z) cos(x).*cos(z).*fIAM(acosd(cos(x).*cos(z)));    % cos(ia)·fIAM(ia)
    lineIAM = @(z) integral(@(x) df(x,z),0,pi/2)/cos(z);
    
    if nargin < 2 || numel(tilt) > 2*(90/DZ)
    % ... pre-calculate an interpolation set, for performance
        zz  = linspace(0,pi/2,ceil(90/DZ));
        interp = true;
    else
        zz = (90 - tilt)*pi/180;
        interp = false;
    end
    
    iam_iso = wedgeIAM(-pi/2,zz); 
    iam_alb = wedgeIAM(zz,pi/2); 
    iam_hb = arrayfun(@(z) lineIAM(z),zz);
    
    if ~interp, varargout = {iam_iso,iam_alb,iam_hb}; return; end
    
    if nargin < 2
        iso = @(T) interp1(zz,iam_iso,(90 - T)*pi/180);
        alb = @(T) interp1(zz,iam_alb,(90 - T)*pi/180);
        hb = @(T) interp1(zz,iam_hb,(90 - T)*pi/180);
    else
        iso = interp1(zz,iam_iso,(90 - tilt)*pi/180);
        alb = interp1(zz,iam_alb,(90 - tilt)*pi/180);
        hb = interp1(zz,iam_hb,(90 - tilt)*pi/180);  
    end
    varargout = {iso,alb,hb};
    
    if nargout > 3, varargout{4} = circumsolarIAM(fIAM,varargin{:}); end
    
    function r = wedgeIAM(a,b)
    % Equivalent to:
    %   r = arrayfun(@(a,b) integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a,b)/(pi/4*(sin(b)-sin(a)));
    %
    % Provisionally solves integration convergence issues for discontinuous functions -- e.g.
    % diffuseIAM('hcpv') -- by lowering tolerance, and capturing warnings.
    %
    % TODO: find robust integration solution for discontinuous functions
    % wedgeIAM = @(a,b) integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a,b)/(pi/4*(sin(b)-sin(a)));

        TOL = 1e-4;
        warning_reseter = naptime();  %#ok<NASGU>
        warning('off','MATLAB:integral2:maxFunEvalsPass');
        warning('error','MATLAB:integral2:maxFunEvalsFail'); %#ok<CTPCT>
        
        [a,b] = compatiblesize(a,b);
        msg = {};
        k = 1./(pi/4*(sin(b)-sin(a)));
        r = NaN(size(k));
        for j = 1:numel(a)
            if a(j) >= b(j), r(j) = 0; continue; end
            try
                r(j) = k(j).*integral2(@(x,z) cos(x).*df(x,z),0,pi/2,a(j),b(j),...
                    'AbsTol',TOL*k(j),'RelTol',TOL);
            catch ERR
                msg{end+1} = getReport(ERR); %#ok<AGROW>
            end
        end
        if ~isempty(msg)
           msg = uniquecell(msg);
           warning(strjoin(msg,newline()));
        end
    end
end

