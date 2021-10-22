function IAM = LI200_IAM(z)
% IAM = LI200_IAM(Z) - cosine response curve for LI-COR sensors, for incidence angle (deg) Z.
%
% NOTE: this is an unofficial custom fit to the curve from the user manual [1], NOT endorsed 
%       by the manufacturer. Values for Z > 85° extrapolated to match IAM = 0 @ 90°.
%
% EXAMPLE: plot(0:90,LI200_IAM(0:90))
%
% [1] LI-COR Terrestrial Radiation Sensors Instruction Manual
% URL: https://www.licor.com/env/support/Light-Sensors/home, retrieved: 19.10.2021.

    p = [0.0237287,0.939568,-2.5646,1.74606];
    q = [-2.60187,1.74562];
    f = @(x) (p(1)*x.^3 + p(2)*x.^2 + p(3)*x + p(4)) ./ (x.^2 + q(1)*x + q(2));
    
    a = 1.25;
    b = 0.032;
    g = @(x) a*exp(-b./cos(x));
    
    x0 = 1.430911022675;
    
    x = min(abs(z),90)*pi/180;
    hi = x > x0;
    IAM = zeros(size(z));
    IAM(hi) = g(x(hi));
    IAM(~hi) = f(x(~hi));
end