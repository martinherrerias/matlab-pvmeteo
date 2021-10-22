function W = point_weights(varargin)
% weights for trapezoidal integration in ND

    W = 1;
    for j = 1:nargin
        x = varargin{j};
        validateattributes(x,{'numeric'},{'real','increasing','vector'});
        if size(x,2) > 1, x = x'; end
        w = movmean([0;diff(x);0],2,'endpoints','discard');
        w = shiftdim(w,1-j);
        W = W.*w;
    end
end