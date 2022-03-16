function TPM = mixed_models(TPM,X)
% TPM = MIXED_MODELS(TPM,K) - Generate model combinations from the output TPM of KNKD_DENSITY_FIT.
%   If TPM contains N 'atomic' models (i.e. rows with non-empty TPM.interpolant), MIXED_MODELS 
%   will make sure to include all relevant* combinations out of the NCHOOSEK(N,K) possibilities. 
%   E.g. if TPM contains models 'a','b','c', MIXED_MODELS(TPM,2) will ensure that '(a)x(b)',
%   '(b)x(c)', and '(a)x(c)' are included in TPM.
%   (*) Not added are combinations with the 'uniform' distribution, or those that contain two or 
%   more static distributions ('constant','separation','copula',..)
%
% TPM = MIXED_MODELS(TPM,TESTS) - Search for, and if missing, add specific model combinations.
%
% See also: KNKD_DENSITY_FIT

    narginchk(2,2);
    if isnumeric(X)
        N = X;
        validateattributes(N,'numeric',{'integer','scalar','>',1});
        TESTS = {};
    else
        TESTS = cellstr(X);
    end

    assert(istable(TPM) && ...
        all(ismember({'conditions','interpolant','texlbl'},TPM.Properties.VariableNames)),'Bad TPM');
    TPM.keys = TPM.Properties.RowNames;
    
    if isempty(TESTS)
        atomic = ~cellfun(@isempty,TPM.interpolant) & ~contains(TPM.keys,'uniform');
        N = min(N,numel(atomic));
        
        static = cellfun(@isempty,TPM.conditions(atomic)); % constant, copula, separation, ...
        
        keys = cellfun(@(x) ['(' x ')'],TPM.keys(atomic),'unif',0);
        TESTS = cell(N-1,1);
        for j = 2:N
            C = nchoosek(1:nnz(atomic),j);
            C(sum(static(C),2) > 1,:) = []; % discard combinations with more than one static dist.
            TESTS{j-1} = arrayfun(@(r) strjoin(keys(C(r,:)),'x'),1:size(C,1),'unif',0);
        end
        TESTS = cat(2,TESTS{:})';
    end
    
    ia = is_synonym(TESTS,TPM.keys);
    TESTS(ia) = [];
    j0 = size(TPM,1);
            
    N_tests = numel(TESTS);
    failed = false(N_tests + j0,1);
    
    warning_disabler = naptime('MATLAB:table:RowsAddedExistingVars'); %#ok<NASGU>

    for j = 1:N_tests
        k = j + j0;
        try
            elements = strsplit(TESTS{j},'x');
            elements = regexprep(elements,'\((.*)\)','$1');

            [~,ie] = ismember(elements,TPM.keys);
            assert(all(ie > 0),'Combination TEST does not match precalculated elements');

            condk = TPM.conditions(ie);
            allcond = unique(cat(2,condk{:}),'stable');

            TPM.keys{k} = TESTS{j};
            % TPM.interpolant{j} = G;
            TPM.conditions{k} = allcond;

            lbl = regexp(TPM.texlbl(ie),'\$(.*) = f(\(.*\))\$','tokens');
            while iscell(lbl{1}), lbl = cat(1,lbl{:}); end
            lbl(:,2) = arrayfun(@(f,a) [f,a{:}],char('e'+(1:size(lbl,1)))',lbl(:,2),'unif',0);
            TPM.texlbl{k} = ['$' lbl{1} ' = ' strjoin(lbl(:,2),' \\cdot ') '$'];

        catch ERR
            failed(k) = true;
            warning('TEST: %s failed with error: %s',TPM.keys{k},ERR.message);
        end        
    end
    TPM.Properties.RowNames = TPM.keys; TPM.keys = [];
    TPM(failed,:) = [];
    
end