classdef Marginal
% An object of class Marginal defines an operation fcn(x,y) over the two-state vector x,y,
% such that m = obj.get(Pxy), where Pxy is a 2D gridded approximation of the PDF of (x,y),
% returns a 1D gridded approximation of the marginal m of fcn(x,y).
%
% TODO: integrate with GriddedPDF class

    properties (GetAccess = public, SetAccess = private)
        grids (1,3) cell
        fcn (1,1) function_handle = @sum
    end
    properties (Hidden = true, GetAccess = public, SetAccess = private)
        W (:,:) double
    end
    properties (Transient = true, GetAccess = public, SetAccess = private)
        gridsize (1,:) double {mustBePositive,mustBeInteger}
        weights (:,1) double {mustBeNonnegative}
    end
    methods
        function obj = Marginal(grids,fcn)
            
            cellfun(@(c) validateattributes(c,'numeric',{'vector','real','increasing'}),grids);

            obj.grids = grids;
            obj.fcn = fcn;
            obj = Marginal.loadobj(obj);
            
            w_xy = point_weights(grids{1:2});
            [xg,yg] = ndgrid(grids{1:2});
            z_qry = fcn(xg(:),yg(:));
            obj.W = interpmatrix(z_qry,grids{3}','extrap','nearest');
            obj.W = obj.W.*double(w_xy(:));
        end
        function Pz = get(obj,Pxy)
            Pz = obj.W'*double(Pxy(:));
            Pz = single(Pz/sum(Pz.*obj.weights));
        end
    end
    methods (Static = true)
        function obj = loadobj(obj)
            obj.gridsize = cellfun(@numel,obj.grids); 
            obj.weights = point_weights(obj.grids{3});
        end
    end
end
