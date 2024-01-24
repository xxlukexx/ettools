function [bins, edges, tall, wide, vecProp, vecRatio] =...
    etBinAOI(in, gaze, binWidthS, piMissing, def, label)

    if ~exist('label', 'var') || isempty(label)
        label = '';
    end

    numAOIs = size(in, 3);
    numSubs = size(in, 2);

    % find edges in gaze time vector
    numBins = round(gaze.Duration / binWidthS);
    edges = 0:binWidthS:binWidthS * (numBins - 1);
    
    % preallocate 
    bins = cell(1, numBins);
    tmp_tall = cell(1, numBins);
    tmp_wide = cell(1, numBins);
    vecProp = nan(numBins, numSubs, numAOIs);
    vecRatio = nan(numBins, numSubs, numAOIs);
    tall = [];
    wide = [];

    if numBins == 0, return, end
    
    % convert bin edges to samples
    edges_samp = arrayfun(@(x) find(gaze.Time >= x, 1), edges);
    s1 = edges_samp;
    s2 = [edges_samp(2:end) - 1, gaze.NumSamples];
    
    % process
    for b = 1:numBins
        
        % sub-segment in variable for this bin
        in_ss = in(s1(b):s2(b), :, :);
        
        % sub-segment gaze
        gaze_ss = gaze.SegmentBySample(s1(b), s2(b));
        
        % sub-segment post-interp missing
        piMissing_ss = piMissing(s1(b):s2(b));
        
        % summarise
        label_ss = sprintf('%s#%d', label, b);
        [bins{b}, tmp_tall{b}, tmp_wide{b}] =...
            etSummariseAOI(in_ss, gaze_ss, def, label_ss, piMissing_ss);
        
        for a = 1:numAOIs
            vecProp(b, :, a) = bins{b}(a).propInAOI;
            vecRatio(b, :, a) = bins{b}(a).ratioInAOI;
        end
        
    end
    
    % flatten tables
    tall = vertcat(tmp_tall{:});
    wide = vertcat(tmp_wide{:});
    
    % extract bin numbers from label
    parts = cellfun(@(x) strsplit(x, '#'), tall.id, 'UniformOutput', false);
    numSepsInOrig = length(strsplit(label, '#'));
    tall.bin = cellfun(@(x) str2double(x{numSepsInOrig + 1}), parts);
    tall.binEdgeS(:, 1) = edges(tall.bin);
    tall = movevars(tall, {'bin', 'binEdgeS'}, 'after', 'id');
    tall.id = repmat({label}, size(tall, 1), 1);
    
    parts = cellfun(@(x) strsplit(x, '#'), wide.id, 'UniformOutput', false);
    wide.bin = cellfun(@(x) str2double(x{numSepsInOrig + 1}), parts);
    wide.binEdgeS(:, 1) = edges(wide.bin);
    wide = movevars(wide, {'bin', 'binEdgeS'}, 'after', 'id');    
    wide.id = repmat({label}, size(wide, 1), 1);

end