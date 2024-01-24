function [smry, looks] = etScoreAOILooks(in, gaze, smry)

    % if no existing summary passed, create an empty struct 
    if ~exist('smry', 'var') || isempty(smry)
        smry = struct;
    end
    
    % only segment gaze data around looks if requested (since it is slow)
    returnLooks = nargout == 2;
    
    % loop through AOIs and extract metrics
    numAOIs = size(in, 3);
    numSubs = size(in, 2);
    
    % loop through each aoi, and each subject
    looks = cell(numAOIs, numSubs);
    for a = 1:numAOIs
        
        for s = 1:numSubs
            
            % find contiguous runs of looking to computer look stats (mean
            % look etc.)
            ct = findcontig2(in(:, s, a));

            % convert from samples to time
            ctt = contig2time(ct, gaze.Time);
            smry(a, s).numLooks = size(ctt, 1);

            % if no continuous runs, make a fake contiguous table with NaNs
            % in it, so that the following stats also return NaN
            if isempty(ctt)
                ctt = [nan, nan, nan];
                looksFound = false;
            else
                looksFound = true;
            end

            % mean look is the mean duration from the ct table
            smry(a, s).meanLook = mean(ctt(:, 3));

            % find the index in the ct table of the peak look
            idx_peakLook = find(ctt(:, 3) == max(ctt(:, 3)), 1);
            smry(a, s).peakLook = ctt(idx_peakLook, 3);

            % since we used find, this will return empty if there were no
            % looks (i.e. the ct table is all NaN). Detect this and
            % overwrite the empty with a NaN
            if isempty(smry(a, s).peakLook)
                smry(a, s).peakLook = nan;
            end

            % do the same for the minimum look
            idx_minLook = find(ctt(:, 3) == min(ctt(:, 3)), 1);
            smry(a, s).minLook = ctt(idx_minLook, 3);
            if isempty(smry(a, s).minLook)
                smry(a, s).minLook = nan;
            end            
            
            % optionally segment data around looks
            if returnLooks && looksFound
                gaze_tmp = gaze.FilterOneSubject(s);
                looks{a, s} = gaze_tmp.SegmentBySample(ct(:, 1), ct(:, 2));
                if ~iscell(looks{a, s}), looks{a, s} = looks(a, s); end
            else 
                looks{a, s} = [];
            end
            
        end
        
    end
        
end