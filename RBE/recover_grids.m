function [grids,var] = recover_grids(TPM,var) 
% recover grid vectors from interpolants

    if nargin < 2, var = uniquecell(cat(2,'kn','kd',TPM.conditions{:}),'stable'); end
    
    M = numel(var);
    grids = cell(1,M);
    for j = 1:size(TPM,1)
        if isempty(TPM.interpolant{j}), continue; end
        [ia,ib] = ismember(var,[{'kn','kd'},TPM.conditions{j}]);
        assert(nnz(ia) == numel(TPM.conditions{j})+2,'Unexpected variable names');

        existing = ~cellfun(@isempty,grids);
        new = ia & ~existing;
        if any(new)
            grids(new) = TPM.interpolant{j}.GridVectors(ib(new));
        end
        old = ia & existing;
        if any(old)
            assert(isequal(grids(old),TPM.interpolant{j}.GridVectors(ib(old))),...
                'Inconsistent GRIDS');
        end
    end 
end