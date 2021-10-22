function [lia,locb] = is_synonym(A,B)
% [lia,locb] = IS_SYNONYM(A,B) - Works like ISMEMBER for cell arrays A, B, except when these
%   can be decomposed into 'atomic' tokens combined with operators 'x' and '+' as used to name
%   transition models in KNKD_DENSITY_FIT. Then conmuted expressions are considered equivalent,
%   e.g. '(a)x(b+c)' == '(c+b)x(a)'.
%
%   Nested parenthesis probably won't work, and there is no distributive propperty, i.e. a+b is
%   considered an inseparable entity (a poor choice of symbols, I know).

    % if nargin == 0, test(); end

    A = cellstr(A);
    B = cellstr(B);
    if isempty(A), lia = []; locb = []; return; end
    if isempty(B), lia = zeros(size(A)); locb = lia; return; end
        
    list = [A(:);B(:)];
    list = strrep(list,' ','');
    
    tokens = regexp(list,'\(([^()]*)\)','tokens','all');
    tokens = cat(2,tokens{:});
    if ~isempty(tokens)
        tokens = cat(2,tokens{:});
        keys = sortedexpr(tokens);
        [keys,~,idx] = uniquecell(keys);
        keys(2,:) = arrayfun(@(j) ['#' num2str(j)],1:numel(keys),'unif',0);
        keys = keys';
        for j = 1:numel(tokens)
            list = strrep(list,['(' tokens{j} ')'],keys{idx(j),2});
        end
        tokens = keys(:,1);
        keys = keys(:,2);
    end
    
	list = sortedexpr(list);
    if ~isempty(tokens)
        for j = 1:numel(tokens)
            list = strrep(list,keys{j},tokens{j});
        end
    end
    
    A(:) = list(1:numel(A));
    B(:) = list(numel(A)+1:end);
    [lia,locb] = ismember(A,B);
end

function X = sortedexpr(X,symbol)

    if nargin < 2, symbol = ''; end

    if iscell(X)
        X = cellfun(@(x) sortedexpr(x,symbol),X,'unif',0); 
        return;
    end

    if isempty(symbol)
        sortplus = @(x) sortedexpr(x,'+'); %#ok<NASGU>
        X = regexprep(X,'([+\w]*)','${sortplus($1)}');
        X = sortedexpr(X,'x');
        return;
    end

    X = strsplit(X,symbol);
    X = sort(X);
    X = strjoin(X,symbol);
end

% function test()
% 
%     A = {'(a)x(b+c)','a+b','d'};
%     B = {'b+a','(c+b)x(a)','e'};
%     [lia,locb] = is_synonym(A,B)
% end