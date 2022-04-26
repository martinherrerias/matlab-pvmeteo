classdef PDFmovmean < matlab.mixin.Copyable
% An object of class PDFmovmean defines a moving average of window N over a single state variable,
% keeping track of the last N distributions passed to the object with obj.push(P), and estimating
% an obj.average assuming either fully-correlated or independent samples.
%
%     obj = PDFmovmean({grid-vector}, N, correlated*);
%     for j = 1:steps, 
%         ...
%         obj.push(P);     % add gridded representations P of a PDF to the stack
%         Q = obj.average; % get the average PDF assuming correlated* or uncorrelated samples
%     end   
%
% TODO: integrate with GriddedPDF class

    properties (GetAccess = public, SetAccess = public)
        grids (1,:) cell
        correlated (1,1) logical = 0
        N (1,1) double {mustBeNonnegative,mustBeInteger}
    end
    properties (Hidden = true, Transient = true, GetAccess = public, SetAccess = private)
        stack (:,:) double
        last (1,1) double {mustBeNonnegative,mustBeInteger}
        marginals
    end
    properties (Transient = true, GetAccess = public, SetAccess = private)
        gridsize (1,:) double {mustBePositive,mustBeInteger}
        weights double {mustBeNonnegative} % result of point_weights(obj.grids)
    end
    properties (Dependent = true)
        average 
    end
    methods
        function obj = PDFmovmean(grids,N,correlated)
            if isnumeric(grids), grids = {grids}; end
            cellfun(@(c) validateattributes(c,'numeric',{'nrows',1,'real','increasing'}),grids);
            obj.grids = grids;
            obj.N = N;
            if nargin > 2, obj.correlated = correlated; end
            obj = PDFmovmean.loadobj(obj);
        end
        function obj = reset(obj,type,P)
        % OBJ = RESET(OBJ,"nan",PDF) - erase stack (set to NaN), then set the first entry to PDF
        % OBJ = RESET(OBJ,"cpy",PDF) - reset stack to N copies of PDF
        
            assert(all(size(P,1:numel(obj.gridsize)) == obj.gridsize),'Non matching P size')
            switch lower(type)
                case 'cpy'
                    obj.stack = repmat(P(:),1,obj.N);
                case 'nan'
                    obj.stack = NaN(numel(P),obj.N);
                    obj.stack(:,1) = P(:);
                otherwise, error("Unknown reset option");
            end
            obj.last = 1;
        end
        function obj = push(obj,P)
        % OBJ = PUSH(OBJ,PDF) - add PDF to the stack (replacing the oldest PDF)
            obj.last = mod(obj.last,obj.N)+1;
            obj.stack(:,obj.last) = P(:);
        end
        function Q = get.average(obj)
            if obj.correlated
                Q = mean(obj.stack,2,'omitnan');
                if numel(obj.gridsize) > 1, Q = reshape(Q,obj.gridsize); end
                Q = Q./sum(Q.*obj.weights);
            else
                assert(numel(obj.grids) == 1,'Uncorrelated average works only for 1D distributions');
                
                % combine X into pairs, then pairs into fourths, etc.
                m = obj.getmarginal(1,1);
                X = obj.stack(:,~all(isnan(obj.stack),1));
                w = ones(size(X,2),1);
                [X,w] = reducepairs(X,w,m);
                
                % ... then combine resulting sets (unequal weights)
                while numel(w) > 1
                    [w,idx] = sort(w); X = X(:,idx);
                    m = obj.getmarginal(w(1),w(2));
                    X(:,1) = m.get(X(:,1).*X(:,2)');
                    w(1) = w(1) + w(2);
                    [X,w] = reducepairs(X,w,m,(1:numel(w))'==2);
                end
                Q = X./sum(X.*obj.weights);
            end
        end
    end
    methods (Access = private, Hidden = true)
        function m = getmarginal(obj,a,b)
            c = gcd(a,b);
            key = sprintf('m_%d_%d',a/c,b/c);
            if isfield(obj.marginals,key), m = obj.marginals.(key);
            else
                m = Marginal(repmat(obj.grids,1,3),@(x,y) a/(a+b)*x + b/(a+b)*y);
                obj.marginals.(key) = m;
            end
        end
    end
    methods (Static = true)
        function obj = loadobj(obj)
            obj.gridsize = cellfun(@numel,obj.grids); 
            obj.weights = point_weights(obj.grids{:});
            obj.stack = NaN(prod(obj.gridsize),obj.N);
            obj.last = 0;
        end
        
        function test()
             GUIfigure('PDFmovmean:test'); clf(); hold on;
             C = PDFmovmean({-6:0.1:6},3, true);
             U = PDFmovmean({-6:0.1:6},3, false);
             for j = 1:10
                 mu = randn(1);
                 sigma = rand(1) + 0.2;
                 P = normpdf(C.grids{1},mu,sigma);
                 C.push(P);
                 U.push(P);
                 
                 if j > 1, delete(h); delete(Hc); delete(Hu); end
                 h = plot(C.grids{1},C.stack,'b-','LineWidth',0.2);
                 Hc = plot(C.grids{1},C.average,'r-','LineWidth',2);
                 Hu = plot(C.grids{1},U.average,'g-','LineWidth',2);
                 pause(1);
             end
        end
    end
end

function [X,w] = reducepairs(X,w,m,used)
% combine X into pairs, then pairs into fourths, etc. until there are no two X with equal w.

    if nargin < 4, used = false(numel(w),1); end
    
    while numel(w) > 1
        didsomething = false;
        for j = 1:numel(w)
           if used(j), continue; end
           k = find(w(j+1:end) == w(j) & ~used(j+1:end),1);
           if ~isempty(k)
               k = k+j;
               X(:,j) = m.get(X(:,j).*X(:,k)');
               w(j) = w(j) + w(k);
               used(k) = true;
               didsomething = true;
           end
        end
        w(used) = [];
        X(:,used) = [];
        used(used) = [];
        if ~didsomething, break; end
    end
end