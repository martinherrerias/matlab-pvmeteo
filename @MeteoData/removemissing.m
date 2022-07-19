function MD = removemissing(MD,fld,E,W,quiet)
% MD = REMOVEMISSING(MD,[FLD,E,W,QUIET]) - remove any fields/columns of MD for which the
%   fraction of flagged points (ignoring MD.missing | MD.dark) exceeds the threshold E 
%   (default MD.options.outlier_error).
%   Issue warnings for fields/columns that exceed threshold W (MD.outputs.outlier_warning).
%
% FLD can be used to restrict the search to certain field(s).
% QUIET = true supresses any warnings or dialogs.

    if nargin < 2 || (isempty(fld) && ~iscell(fld)), fld = fieldnames(MD); 
    else
        [fld,ix] = parselist(fld,fieldnames(MD));
    end

    if isempty(MD.missing) || isempty(MD.dark)
        filter = true(MD.Nt,1);
        if nargin < 3 || isempty(E), E = 1; end
        if nargin < 4 || isempty(W), W = 1; end
    else
        filter = ~MD.missing & ~MD.dark;
        if nargin < 3 || isempty(E), E = MD.options.outlier_error; end
        if nargin < 4 || isempty(W), W = MD.options.outlier_warning; end
    end
    if nargin < 5, quiet = false; end

    n = nnz(filter);
    if n == 0
        warning('All steps flagged as missing/dark');
        return; 
    end

    for j = 1:numel(fld)
%         available = repmat(filter,1,size(MD.data.(fld{j})));
%         available(filter,:) = ~isnan(MD.data.(fld{j})(filter,:));
%         na = sum(available,1);
%         
        bad_fraction = sum(MD.flags.data.(fld{j})(filter,:) > 0,1)/n;
        src = MD.getsourceof(fld{j});

        fishy = bad_fraction >= W & bad_fraction < E;
        if any(fishy) && ~quiet
            msg = flagsummary(MD.flags,fld{j},filter);
            if ~all(fishy) 
                msg = strsplit(msg,newline());
                msg = strjoin(msg(fishy),newline());
                msg = sprintf('You might want to check %s.%s (%s):\n%s',...
                              fld{j},shortliststr(find(fishy)),...
                              shortliststr(src(fishy)),msg);
            else
                msg = sprintf('You might want to check field %s:\n%s',fld{j},msg);
            end
            warning('meteodata:outlier_warning','%s',msg);
        end

        bad = bad_fraction >= E;
        if ~any(bad), continue; end
        if ~quiet
            msg = MD.flags.flagsummary(fld{j},filter);
            if ~all(bad)
                msg = strsplit(msg,newline());
                msg = strjoin(msg(bad),newline());
                msg = sprintf('Flagged %s.%s (%s) will be removed from meteo-data:\n%s',...
                              fld{j},shortliststr(find(bad)),shortliststr(src(bad)),msg);
            else
                msg = sprintf('Flagged field %s will be removed from meteo-data:\n%s',...
                              fld{j},msg);
            end

            if runningfromUI()
                switch optquestdlg([msg ', proceed?'],'All-flagged','Yes','What? No!','Yes')
                    case 'Yes'
                    case 'What? No!', continue; 
                    otherwise, error('User cancelled diaglog');
                end
            else
                warning('cmd:nanfield',msg);
            end
        end

        if all(bad)
            MD.data.(fld{j}) = [];
            MD.data.Properties.CustomProperties.source{ix(j)} = [];
            MD.flags.data.(fld{j}) = [];
        elseif any(bad)  
            MD.flags.flagsummary(fld{j},filter)
            MD.data.(fld{j})(:,bad) = [];
            MD.data.Properties.CustomProperties.source{ix(j)}(bad) = [];
            MD.flags.data.(fld{j})(:,bad) = [];
        end
    end 
end