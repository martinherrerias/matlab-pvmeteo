function mcclear_healthcheck(FULL)
% mcclear_healthcheck(FULL) - read headers (or full files) to check for errors. 
%   Standardize file names for faster searching. 
%
% TODO: look for overlaps and delete redundant files

    if nargin < 1, FULL = false; end
    EXT = '.csv';
    PREFIX = 'mcclear';

    files = pickfile(['*' EXT],Inf);

    valid = true(size(files));
    stdnames = cell(size(files));
    for j = 1:numel(files)
        try 
            F = read_CAMS(files{j},'data',FULL);
            stdnames{j} = meteofilename(F,F.date_begin,F.date_end,F.summarization,'prefix',PREFIX);
        catch ERR
            valid(j) = false;
        end
    end
    if ~all(valid)
        warning('Failed to read %s. Renaming to *.bak',nthings(nnz(~valid),'file'));
        cellfun(@movefile,files(~valid),regexprep(files(~valid),['(\' EXT ')$'],'.bak'));
    end

    [paths,files] = cellfun(@fileparts,files,'unif',0);
    stdz = strcmp(files,stdnames);
    stdz(~valid) = true;

    if ~all(stdz)
    % offer to rename to standardized names, for faster evaluation next time
        msg = sprintf('%s not named to standard. Rename?',nthings(nnz(~stdz),'file'));
        switch optquestdlg(msg,'healthckeck:rename','yes','no','no')
        case 'yes'
            files = cellfun(@(p,f) fullfile(p,[f EXT]),paths,files,'unif',0);
            stdnames = cellfun(@(p,f) fullfile(p,[f EXT]),paths,stdnames,'unif',0);
            cellfun(@movefile,files(~stdz),stdnames(~stdz)); 
        end
    end

end









