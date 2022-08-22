classdef GriddedPDF < matlab.mixin.Copyable
% GRIDDEDPDF - class to represent N-dimensional Probability Distribution Functions, sampled over
% irregular rectangular grids. Each sample (grid point) is weighted according to the average
% distance to its neighbors, so that the sum of weighted probabilities sums to 1 (trapezoidal
% integration).
%
% TODO: rewrite RBE functions in terms of this class.

    properties (GetAccess = public, SetAccess = public, AbortSet)
        P double % {mustBeNonnegative}
    end
    properties (GetAccess = public, SetAccess = private)
        grids (1,:) cell
        lastnorm (1,1) double {mustBeNonnegative}
    end
    properties (Transient = true, GetAccess = public, SetAccess = private)
        gridsize (1,:) double {mustBePositive,mustBeInteger} = []
        weightvec (1,:) cell
        dims (1,1) double {mustBePositive,mustBeInteger} = 1
        weights double {mustBeNonnegative}
        mean
    end
    properties (Transient = true, Hidden = true) % Calculated (once) on demand
        Dxx
        wP double
        marginals
    end
    properties (Dependent = true)
        mode
        median
        cov
        ndgrids cell
    end
    methods
        function obj = GriddedPDF(grids,P)
            if isnumeric(grids), grids = {grids}; end
            cellfun(@(c) validateattributes(c,'numeric',{'nrows',1,'real','increasing'}),grids);
            obj.grids = grids;
            obj = GriddedPDF.loadobj(obj);
            obj.P = P;
        end
        function set.P(obj,P)
            if isempty(obj.gridsize) || isempty(obj.weights) || isempty(obj.weights) %#ok<MCSUP>
                obj = GriddedPDF.loadobj(obj); % solves initialization issue (<MCSUP> warning)
                justloaded = true;
            else
                justloaded = false;
            end
            P = compatiblesize(P,obj.weights);    %#ok<MCSUP>
            w = sum(P.*obj.weights,'all');        %#ok<MCSUP>
            if ~justloaded, obj.lastnorm = w; end %#ok<MCSUP>
            if w > 0
                obj.P = P/w;
            else
                if isempty(obj.P), obj.P = zeros(obj.gridsize); end %#ok<MCSUP>
                obj.P(:) = 1/sum(obj.weights,'all'); %#ok<MCSUP>
            end
            obj.reset();
        end
        function reset(obj)
            obj.wP = [];
            obj.mean = [];
            obj.marginals = [];
        end
        function g = get.ndgrids(obj), [g{1:obj.dims}] = ndgrid(obj.grids{:}); end
        function W = get.wP(obj)
            if isempty(obj.wP)
                W = obj.weights.*obj.P; 
                obj.wP = W;
            else
                W = obj.wP;
            end 
        end
        function D = get.Dxx(obj)
            assert(maxarraysize()/prod(obj.gridsize) > 2,'Out of memory')
            if isempty(obj.Dxx)
                g = cellfun(@(x) reshape(x,[],1),obj.ndgrids,'unif',0);
                Xg = cat(2,g{:});
                obj.Dxx = squareform(pdist(Xg,'euclidean'));
            end
            D = obj.Dxx; 
        end
        function S = stats(obj)
            S = struct('mean',obj.mean,'mode',obj.mode,'median',obj.median,'cov',obj.cov);
        end
        function P = marginal(obj,dim)
            assert(obj.dims > 1 && dim <= obj.dims);
            if isempty(obj.marginals), obj.marginals = cell(1,obj.dims); end
            if isempty(obj.marginals{dim})
                w = obj.weightvec{dim};
                p = shiftdim(sum(obj.wP,[1:dim-1,dim+1:obj.dims]),dim-1)./w;
                obj.marginals{dim} = GriddedPDF(obj.grids(dim),p);
            end
            P = obj.marginals{dim};
        end
        function [X,C] = cdf(obj)
        % [X,C] = CDF(OBJ) - return points X and values C that form a piece-wise-linear approx.
        %   to the Cumulative Distribution Function for OBJ. The points X are midpoints of the
        %   original grid, plus the first and last points.
        
            X = cell(1,obj.dims);
            C = obj.wP;
            for j = 1:obj.dims
                midpts = interpn(obj.grids{j},1);
                X{j} = midpts([1,2:2:end,end]);
                s = obj.gridsize;
                s(j) = 1;
                C = cat(j,zeros(s),cumsum(C,j));
            end
            if numel(X) == 1, X = X{1}; end
        end
        function m = get.median(obj), m = quantile(obj,0.5); end
       
        function q = quantile(obj,p,dim)
            if nargin < 3, dim = 1:obj.dims; end
            q = zeros(numel(p),numel(dim));
            for j = dim(:)'
               [x,c] = marginal(obj,j).cdf();
               [c,idx] = unique(c,'sorted');
               q(:,j) = interp1(c,x(idx),p);
            end
        end
        function m = get.mean(obj)
            if isempty(obj.mean)
                m = zeros(1,obj.dims);
                for k = 1:obj.dims
                    d = [1:k-1,k+1:obj.dims];
                    if isempty(d), w = obj.wP; else, w = sum(obj.wP,d); end
                    m(k) = sum(shiftdim(obj.grids{k},2-k).*w,k);
                end
                obj.mean = m;
            end
            m = obj.mean;
        end
        function C = get.cov(obj)
            m = obj.mean;
            C = zeros(obj.dims);
            for i = 1:obj.dims
                for j = i:obj.dims
                    d = [1:i-1,i+1:j-1,j+1:obj.dims];
                    if isempty(d), w = obj.wP; else, w = sum(obj.wP,d); end
                    C(i,j) = sum(shiftdim(obj.grids{i}-m(i),2-i).*...
                                 shiftdim(obj.grids{j}-m(j),2-j).*w,'all');
                    if j > i, C(j,i) = C(i,j); end
                end
            end
        end
        function m = get.mode(obj)
            [~,k] = max(obj.P,[],'all','linear');
            [i{1:obj.dims}] = ind2sub(obj.gridsize,k);
            m = cellfun(@(c,i) c(i),obj.grids,i);
        end
        function px = interp(obj,varargin)
            if numel(varargin) == 1 && size(varargin,2) == obj.dims
                validateattributes(varargin,'numeric',{'real','size',[NaN,obj.dims]});
                varargin = num2cell(varargin,1);
            else
                s = size(varargin{1});
                cellfun(@(c) validateattributes(c,'numeric',{'real','size',s}),varargin);
            end
            px = interpn(obj.grids{:},obj.P,varargin{:});
        end
        function p = PIT(obj,x,dim)
        % Probability Integral Transform
            if nargin < 3 && obj.dims == 1, dim = 1; else, narginchk(3,3); end
            [g,c] = marginal(obj,dim).cdf();
            p = interp1(g,c,x);
        end
        function [b,px] = BOT(obj,varargin)
        % Box Ordinate Transform
            px = interp(obj,varargin{:});
            b = 1-sum(obj.wP(obj.P <= px),'all');
        end
        function e = energyscore(obj,varargin)
        % Energy Score
            if numel(varargin) == 1 && size(varargin,2) == obj.dims
                validateattributes(varargin,'numeric',{'real','size',[1,obj.dims]});
                x = varargin;
            else
                cellfun(@(c) validateattributes(c,'numeric',{'real','scalar'}),varargin);
                x = cat(2,varargin{:});
            end
            g = cellfun(@(x) reshape(x,[],1),obj.ndgrids,'unif',0);
            Xg = cat(2,g{:});
            W = reshape(obj.wP,numel(obj.weights),[])';
            Dxxi = single(vecnorm(Xg - x,2,2));
            e = W*Dxxi - 0.5*sum(W.*(W*obj.Dxx),2);
        end
        function B = resample(A,n,varargin)
        % B = RESAMPLE(A,N,..) - refine grids by N, interpolate values
            B.grids = cellfun(@(x) interpn(x,n,varargin{:}),A.grids,'unif',0);
            [g{1:A.dims}] = ndgrid(B.grids{:});
            B.P = interpn(A.grids{:},A.P,g{:},varargin{:});
            B = GriddedPDF.loadobj(B);
        end
    end
    methods (Static = true)
        function obj = loadobj(obj)
            if ~isa(obj,'GriddedPDF')
                obj = GriddedPDF(obj.grids,obj.P);
                return;
            end
            obj.gridsize = cellfun(@numel,obj.grids); 
            obj.weightvec = cellfun(@point_weights,obj.grids,'unif',0);
            obj.dims = numel(obj.gridsize);
            obj.weights = obj.weightvec{1};
            for j = 2:obj.dims
                obj.weights = obj.weights.*shiftdim(obj.weightvec{j},1-j);
            end
        end
        function varargout = test()
            
            g = {linspace(-6,6,50),linspace(-3,3,40)};
            MU = [1,0.5];
            SIGMA = [1,0.3;0.3,0.5];
            
            [x,y] = ndgrid(g{:});
            P = mvnpdf([x(:),y(:)],MU,SIGMA);
            P = reshape(P,size(x));
            obj = GriddedPDF(g,P);
            
            GUIfigure('GriddedPDF:test'); clf();
            subplot(3,3,[1,2,4,5]); hold on;
                contour(obj.grids{:},obj.P',10);
                pts = [obj.mean; obj.median; obj.mode];
                plot(pts(:,1),pts(:,2),'o');
                lims = axis();
            subplot(3,3,[3,6]);
                my = marginal(obj,2);
                plot(my.P,my.grids{:});
                ylim(lims(3:4));
            subplot(3,3,7:8);
                mx = marginal(obj,1);
                plot(mx.grids{:},mx.P);
                xlim(lims(1:2));
                
            if nargout > 0, varargout{1} = obj; end
        end
    end
end