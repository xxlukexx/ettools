function [gaze, in, res] = etMakeShamData(numSubs, duration, fs, file_aoi, aoi_def)
    
% setup

    numSamps = duration * fs;

% find AOI centroids

    % load AOI
    aoi = imread(file_aoi);
    numAOIs = size(aoi_def, 1);
    
    % binarise image into separate logical mask for each AOI
    bin = etBinariseAOIMask(aoi, aoi_def);

    % find AOI centroids
    [cx, cy] = etFindAOICentroid(bin);
    
    
% plan quantity of data

    propAOI = rand(numSubs, numAOIs);
    
    % noramlise probabilities so they add up to 1 for each AOI
    for r = 1:size(propAOI, 1)
        propAOI(r, :) = propAOI(r, :) ./ sum(propAOI(r, :));
    end
    
    % convert prop to samples
    sampsAOI = floor(propAOI * numSamps);
    
% insert AOI looking

    lx = nan(numSamps, numSubs);
    ly = nan(numSamps, numSubs);
    rx = nan(numSamps, numSubs);
    ry = nan(numSamps, numSubs);
    inAnyAOI = zeros(numSamps, numSubs);
    
    for s = 1:numSubs
        
        s1 = 1;
        for a = 1:numAOIs
            
            % randomise AOI order
            ord = randperm(numAOIs);
            
            s2 = s1 + sampsAOI(s, a) - 1;
            lx(s1:s2, s) = repmat(cx(ord(a)), sampsAOI(s, a), 1);
            ly(s1:s2, s) = repmat(cy(ord(a)), sampsAOI(s, a), 1);
            rx(s1:s2, s) = repmat(cx(ord(a)), sampsAOI(s, a), 1);
            ry(s1:s2, s) = repmat(cy(ord(a)), sampsAOI(s, a), 1);
            inAnyAOI(s1:s2, s) = ord(a);
            
            s1 = s2 + 1;
            
        end
        
    end
    
% plan missing data

    propMissing = 0.1 + (rand(numSubs, 1) * 0.5);
    sampsMissing = round(propMissing * numSamps);
    
    % define distribution of gap length
    gap_mean = 0.30 * fs;
    gap_sd = 1 * fs;
    
% insert missing data. For each subject, loop, inserting chunks of missing
% data at random samples until the total number of missing samples has been
% reached

    missing = false(numSamps, numSubs);
    for s = 1:numSubs
        
        currentMissing = sum(missing(:, s));
        missingNeeded = sampsMissing(s) - currentMissing;
        while missingNeeded > 0
        
            % calculate amount of time missing in seconds
            gapLengthSamps = 1 + round(abs(normrnd(gap_mean, gap_sd)));
            s1 = randi(numSamps);
            s2 = s1 + gapLengthSamps - 1;
            if s2 > numSamps, s2 = numSamps; end
            if s2 - s1 > missingNeeded
                s2 = s1 + missingNeeded - 1;
            end
            
            missing(s1:s2, s) = true;
            currentMissing = sum(missing(:, s));
            missingNeeded = sampsMissing(s) - currentMissing;
        
        end
            
    end

% plan absent data. Everyone has data at the start, but a random 0-20% at
% the end is missing

    propAbsent = rand(numSubs, 1) * .2;
    sampsAbsent = round(propAbsent * numSamps);
    
    % create absent matrix
    absent = false(size(missing));
    for s = 1:numSubs
        s2 = numSamps;
        s1 = s2 - sampsAbsent(s);
        absent(s1:s2, s) = true;
    end

    % make absent gaze samples NaN
    lx(absent) = nan;
    ly(absent) = nan;
    rx(absent) = nan;
    ry(absent) = nan;
    inAnyAOI(absent) = 0;
    
    % absent samples cannot be missing
    missing(absent) = false;
    
    % insert missing data
    lx(missing) = nan;
    ly(missing) = nan;
    rx(missing) = nan;
    ry(missing) = nan;
    inAnyAOI(missing) = 0;
    
% make time vector
    
    sps = 1 / fs;
    t(:, 1) = sps:sps:duration;
    
% make gaze object
    
    gaze = etGazeDataBino;
    gaze.Import(lx, ly, rx, ry, t, missing, missing, absent);
    
% figure out final numbers

    res.propValid = sum(~missing & ~absent, 1) ./ sum(~absent, 1);
    res.sampsAOI = nan(numSubs, numAOIs);
    res.propAOI = nan(numSubs, numAOIs);
    in = false(numSamps, numSubs, numAOIs);
    for a = 1:numAOIs 
        res.sampsAOI(:, a) = sum(inAnyAOI == a, 1);
        res.propAOI(:, a) = res.sampsAOI(:, a) ./ (sum(~missing & ~absent, 1)');
        in(:, :, a) = inAnyAOI == a;
    end












% % setup 
% 
%     % calculate number of samples
%     numSamps = duration * fs;
%     m = nan(numSamps, numSubs, 2);      % third dim is x, y
%     
%     % remove multi-colour AOIs from def
%     isMulti = cellfun(@(x) length(x) ~= 1, aoi_def(:, 2));
%     aoi_def(isMulti, :) = [];
%     
%     % extract colours from aoi def
%     tmp = cellfun(@(x) x{1}, aoi_def(:, 2), 'UniformOutput', false);
%     def_col = vertcat(tmp{:});
%     
% % construct missing data. Onset of each period of missing data is random
% % through the data. Duration of each period if sampled from a probability
% % distribution from 1 sample to half of the total duration. 
% 
%     % number of gaps is random for each subject
%     numGaps = round((0.3 * rand([numSubs, 1])) * duration);
%     
%     % random onset of each gap in data
%     gap_onset = arrayfun(@(x) 1 + round(rand(x, 1) * (numSamps - 1)), numGaps,...
%         'UniformOutput', false);
%     
%     % define distribution of gap length
%     gap_mean = 0.30 * fs;
%     gap_sd = 1 * fs;
%     
%     % sample from normal distribution to find gap lengths
%     gap_length = arrayfun(@(x) 1 + round(abs(normrnd(gap_mean, gap_sd, [x, 1]))),...
%         numGaps, 'UniformOutput', false);
%     
%     % make missing matrix from gap onset and lengths
%     missing = false(numSamps, numSubs);
%     for s = 1:numSubs
%         for g = 1:numGaps(s)
%             
%             s1 = gap_onset{s}(g);
%             s2 = s1 + gap_length{s}(g);
%             if s2 > numSamps, s2 = numSamps; end
%             missing(s1:s2, s) = true;
%             
%          end        
%     end
%     
% % calculate fixation duration per subject, and then allocated number of
% % fixations to generate
% 
%     % fixation duration for each subject
%     fixDur = abs(normrnd(0.300, 0.150, [numSubs, 1]));
%     fixDur_samps = round(fixDur * fs);
%     numFix = round(duration ./ fixDur);
%     
% % for each subject, generate numFix fixations at random locations. If this
% % location is over an AOI then keep it, otherwise generate another
% 
%     % load AOI
%     aoi = imread(file_aoi);
%     numAOIs = size(aoi_def, 1);
%     w = size(aoi, 2);
%     h = size(aoi, 1);
%     
%     % prepare in vector
%     in = false(numSamps, numSubs, numAOIs);
% 
%     for s = 1:numSubs
%         
%         % start at first sample
%         s1 = 1;
%         while s1 <= numSamps
%             
%             % every 1 in 10 fixations is non-AOI looking. Otherwise, loop
%             % until the fixation is over an AOI
%             fixHappy = false;
%             overAOI = rand < .1;
%             while ~fixHappy
% 
%                 % random fixation location
%                 x = rand;
%                 y = rand;
% 
%                 % convert to pixels
%                 x_px = 1 + round(x * (w - 1));
%                 y_px = 1 + round(y * (h - 1));
% 
%                 % get colour at location
%                 col = shiftdim(aoi(y_px, x_px, :), 2)';
% 
%                 % compare to AOI def
%                 col_dis = abs(col - def_col);
%                 col_match = all(col_dis < 10, 2);
%                 idx_match = find(col_match, 1);
%                 
%                 % decide on whether this has to be over an AOI or not
%                 fixHappy = overAOI || ~isempty(idx_match);
%             
%             end
%             
%             % place fixation in data
%             s2 = s1 + fixDur_samps(s);
%             if s2 > numSamps, s2 = numSamps; end
%             m(s1:s2, s, 1) = repmat(x, s2 - s1 + 1, 1);
%             m(s1:s2, s, 2) = repmat(y, s2 - s1 + 1, 1);
%             
%             % put AOI stats in in vector
%             if overAOI
%                 in(s1:s2, idx_match) = true;
%             end
%             
%             s1 = s2 + 1;
%             
% %             tmpMiss = false(size(missing));
% %             tmpMiss(s1:s2, s) = missing(s1:s2, s);
% %             m(tmpMiss, :) = nan;
% %             
%         end
%         
%     end
%     
% % replace missing data
% 
%     tmp_x = m(:, :, 1);
%     tmp_y = m(:, :, 2);
%     tmp_x(missing) = nan;
%     tmp_y(missing) = nan;
%     m(:, :, 1) = tmp_x;
%     m(:, :, 2) = tmp_y;
%     
% % convert data matrix to main and time buffers
% 
%     % make vector of timestamps in µsecs
%     usPerSamp = round(1e6 / fs);
%     tbv(:, 1) = 0:usPerSamp:(numSamps - 1) * usPerSamp;
% 
%     % make timebuffers
%     tb = arrayfun(@(x) [tbv, zeros(numSamps, 1)], 1:numSubs, 'UniformOutput', false);
%     
%     % make empty main buffers
%     mb = arrayfun(@(x) nan(numSamps, 26), 1:numSubs, 'UniformOutput', false);
%     
%     for s = 1:numSubs
%         mb{s}(:, [7, 20]) = repmat(m(:, s, 1), 1, 2);
%         mb{s}(:, [8, 21]) = repmat(m(:, s, 2), 1, 2);
%         mb{s}(missing(:, s), [13, 26]) = repmat(4, sum(missing(:, s)), 2);
%         mb{s}(~missing(:, s), [13, 26]) = zeros(sum(~missing(:, s)), 2);
% 
%     end
   
end