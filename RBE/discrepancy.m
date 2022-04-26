function d = discrepancy(x,a,b)
% D = discrepancy(X,[A,B]) - Returns the Kolmogorov-Smirnov distance between the empirical CDF 
%   of X and the uniform distribution in [A,B]

    x = x(:);
    if nargin < 3 || isempty(a) || isempty(b), [a,b] = bounds(x); end
    
    x = sort(x(isfinite(x)));
    n = numel(x);

    if ~(b > a), d = 1; return; end
    if a ~= 0 || b ~= 1, x = (x-a)/(b-a); end

    d = max(abs(x - (1:n)'/n));
end