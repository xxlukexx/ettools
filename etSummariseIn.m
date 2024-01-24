function [smry, tab_tall, tab_wide, looks] = etSummariseIn(type, in, gaze,...
    def, id, postInterpMissing, dist)
% This function is called by etSummariseAOI and etSummariseSalience. They
% pass all inputs along with 'type' which is either 'aoi' or 'salience'. 

    isAOI = strcmpi(type, 'aoi');
    isSal = strcmpi(type, 'salience');
    if ~isAOI && ~isSal
        error('Type must be either ''aoi'' or ''salience''')
    end

% setup and data/arg check

    % progress reporting 
    startTime = cputime;
    wb = [];
    
    % loop through AOIs and extract metrics
    numAOIs = size(in, 3);
    numSubs = size(in, 2);
    smry = cell(numAOIs, numSubs);
    looks = cell(numAOIs, numSubs);
    
    % if no AOI def passed, make generic AOI labels
    if exist('def', 'var') && ~isempty(def)
        % check that the number of AOIs in def match the results
        if size(def, 1) ~= numAOIs 
            error('Number of AOIs in AOI definition (%d) does not match the number of AOIs in score matrix (%d).',...
                size(def, 1), numAOIs)
        end
        % extract labels
        lab = def(:, 1);
        
    else
        % generate labels in the form of AOI_01, etc.
        lab = arrayfun(@(x) sprintf('%s_%01d', upper(type), x), 1:numAOIs,...
            'UniformOutput', false)';
        
    end
    
    % if no IDs passed, generate numerical index in the form of 0001 etc.
    if ~exist('id', 'var') || isempty(id)
        id = arrayfun(@(x) sprintf('%04d', x), 1:numSubs,...
            'UniformOutput', false)'; 
    end
    
    % if ID has been passed as a single char, place it in a cell array
    if ischar(id)
        id = {id};
    end
    
    % check that number of subs matches size of score matrix
    if length(id) ~= numSubs
        error('Number of IDs (%d) does not match number of subjects in score matrix (%d).',...
            length(id), numSubs)
    end
    
    % if post-interp prop valid has not been passed, use prop valid from
    % the gaze data
    if ~exist('postInterpMissing', 'var') || isempty(postInterpMissing)
        propVal = gaze.PropValid;
        missing = gaze.Missing;
    else
        propVal = sum(~postInterpMissing & ~gaze.Absent, 1) ./ sum(~gaze.Absent, 1);
        missing = postInterpMissing;
    end
    
    % if no distance (to AOI centroid) is passed, set to NaN so output is
    % also NaN
    if ~exist('dist', 'var') || isempty(dist)
        dist = nan(size(in));
    end
    
    % if 'looks' is asked for, we need to extract segments of gaze data
    % around each look and return these. This is done in etScoreAOILooks
    % and is slow, so only do it if asked
    returnLooks = nargout == 4;
    
% loop through subjects, then AOIs, and summarise each
    
    % empty tall and wide tables
    tab_tmp_tall = cell(numSubs, 1);
    tab_tmp_wide = cell(numSubs, 1);
    
    for s = 1:numSubs
        
        % show waitbar if more than 2s has elapsed since the function was
        % called
        if cputime - startTime > 2 && mod(s, 50) == 0
            waitMsg = sprintf('Summarising %s [dataset %d of %d]', type, s, numSubs);
            if isempty(wb)
                wb = waitbar(s / numSubs, waitMsg);
            else
                wb = waitbar(s / numSubs, wb, waitMsg);
            end
        end
        
    % calculate all per-AOI stats
    
        % number of valid (not missing and not absent) samples per subject
        subNumValid = sum(~missing(:, s) & ~gaze.Absent(:, s), 1);
        
        % get gaze for just this subject
        gaze_sub = gaze.FilterOneSubject(s);
    
        % loop through each AOI/salience map
        for a = 1:numAOIs
            
            if isAOI && returnLooks
                [smry{a, s}, looks(a, s)] = etScoreAOILooks(in(:, s, a), gaze_sub, smry{a, s});
            elseif isAOI && ~returnLooks
                smry{a, s} = etScoreAOILooks(in(:, s, a), gaze_sub, smry{a, s});
            end
        
            smry{a, s}.id = id{s};
            smry{a, s}.propVal = propVal(s);
            smry{a, s}.totalScreenWatchTime = propVal(s) * gaze.Duration;
            smry{a, s}.aoi = lab{a};
            
            if isAOI
                
                smry{a, s}.triggered = any(in(:, s, a));
                smry{a, s}.samplesInAOI = sum(in(:, s, a));
                smry{a, s}.propInAOI = sum(in(:, s, a)) / subNumValid;
                timeDelta = [0; diff(gaze.Time)];
                smry{a, s}.timeInAOI = sum(timeDelta(in(:, s, a)));
                smry{a, s}.firstSamp = find(in(:, s, a), 1);
                smry{a, s}.firstTimeS = gaze.Time(smry{a, s}.firstSamp);
                if isempty(smry{a, s}.firstTimeS), smry{a, s}.firstTimeS = inf; end
                if isempty(smry{a, s}.firstSamp), smry{a, s}.firstSamp = inf; end
                smry{a, s}.meanDistance = nanmean(dist(:, s, a));

            elseif isSal
                
                smry{a, s}.salience_mean = nanmean(in(:, s, a));
                smry{a, s}.salience_sd = nanstd(in(:, s, a));

            end
            
            % store inAOI vector
            smry{a, s}.in = in(:, s, a);
            
            % store time vector
            smry{a, s}.time = gaze.Time;
            
            % store distance vector
            if ~isempty(dist)
                smry{a, s}.distance = dist(:, s, a);
            else
                smry{a, s}.distance = [];
            end

        end
        
    % loop again to calculate relative stats (e.g. ratio of looking
    % time between AOIs, first look etc.)
    
        % cat summary across AOIs to create subject-summary
        smry_sub = vertcat(smry{:, s});
    
        if isAOI
    
            % get first look time in seconds for each AOI
            firstLook_tmp = [smry_sub.firstTimeS];

            % find index of AOI which was looked at first
            idx_firstLook = firstLook_tmp == min(firstLook_tmp);

            % get total looking time to all AOIs for this subject
            timeTot = sum([smry_sub.timeInAOI]);

            % calculate ratio in AOI and label of first looked-at AOI
            for a = 1:numAOIs
                smry{a, s}.ratioInAOI = smry{a, s}.timeInAOI / timeTot;
                smry{a, s}.firstLook = idx_firstLook(a);
            end
            
        end

        % make table of results
        smry_tmp = rmfield(smry_sub, {'in', 'time', 'distance'});
        tab_tmp_tall{s} = struct2table(smry_tmp, 'AsArray', true);

        % calculate duration of non-absent data
        if all(gaze.Absent(:, s))
            dataDur = 0;
        else
            dataDur = max(gaze.Time(~gaze.Absent(:, s)));
        end
        tab_tmp_tall{s}.dataDuration =...
            repmat(dataDur, size(tab_tmp_tall{s}, 1), 1);
            
    % build a wide table for this subject
    
        % prepare the table
        tab_tmp_wide{s} = table;
        smry_tmp_c = cell(1, numAOIs);
        fnames = cell(numAOIs, 1);
        
        for a = 1:numAOIs
            
            smry_tmp = rmfield(smry{a, s}, {'in', 'time', 'distance'});
            aoiName = smry_tmp.aoi;
            
            % store AOI name, then remove unwanted fields
            smry_tmp =...
                rmfieldIfPresent(smry_tmp, {'id', 'aoi', 'firstLook',...
                    'propVal', 'totalScreenWatchTime'});        

            % rename fields to refer to the current AOI
            fnames{a} = fieldnames(smry_tmp);
            fnames{a} = cellfun(@(x) sprintf('%s_%s', aoiName, x),...
                fnames{a}, 'uniform', false);
            fnames{a} = fixTableVariableNames(fnames{a})';
            smry_tmp_c{a} = struct2cell(smry_tmp)';

        end
        
        % put table together, column-wise
        smry_tmp_c_all = horzcat(smry_tmp_c{:});
        fnames_all = horzcat(fnames{:});
        tab_tmp_wide{s} = cell2table(smry_tmp_c_all, 'VariableNames', fnames_all);
        
        % add ID to temporary table
        tab_tmp_wide{s}.id = id(s);
        
        % store prop valid (this doesn't change per AOI, so we just
        % store it once in the wide table)
        tab_tmp_wide{s}.propVal = propVal(s);
            
        % compute first look
        if isAOI
            tab_tmp_wide{s}.firstLook{1} = smry_sub(idx_firstLook).aoi;
        end
        
        % store data duration
        tab_tmp_wide{s}.dataDuration{1} = dataDur;
        
        % store total watch time
        tab_tmp_wide{s}.totalScreenWatchTime{1} =...
            smry_sub(1).totalScreenWatchTime;
                
    end    
        
    % cat all temp (per-subject) tables into one big table
    tab_tall = vertcat(tab_tmp_tall{:});
    tab_wide = vertcat(tab_tmp_wide{:});
    
    % reorder wide table columns so that ID, propVal and firstLook are at
    % the start
    tab_wide = [tab_wide(:, end - 3:end), tab_wide(:, 1:end - 4)];

    % convert IDs to cell arrays of strings - note this may error with
    % numeric IDS, if it does, add type checking for id column of table (or
    % force even numeric IDs to be strings)
    tab_tall.id = cellstr(tab_tall.id);
    tab_wide.id = cellstr(tab_wide.id);
    
    % ensure dataDuration column is not a cell array
    if iscell(tab_wide.dataDuration)
        tab_wide.dataDuration = cell2mat(tab_wide.dataDuration);
    end
    if iscell(tab_tall.dataDuration)
        tab_tall.dataDuration = cell2mat(tab_tall.dataDuration);
    end
    if iscell(tab_wide.totalScreenWatchTime)
        tab_wide.totalScreenWatchTime = cell2mat(tab_wide.totalScreenWatchTime);
    end
    if iscell(tab_tall.totalScreenWatchTime)
        tab_tall.totalScreenWatchTime = cell2mat(tab_tall.totalScreenWatchTime);
    end    
   
    if ~isempty(wb), delete(wb), end
    
    % convert cell array of summaries to a struct array
    smry = cell2mat(smry);
    
% profile off
% profile viewer
% pause
    
end