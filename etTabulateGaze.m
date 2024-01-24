function gaze = etTabulateGaze(gazeArray)
% takes and array of teGazeData object and tabulates them into one teGaze
% object. Assume sampling rates are equivalent, resample before using this
% function is not. If durations are not equivalent between gaze objects,
% the table is build on the longest duration, and any shorter gaze objects
% have extraneous gaze samples marked as absent

    numGaze = length(gazeArray);

% find lengths

    lens = [gazeArray.NumSamples];
    idx_longest = find(lens == max(lens), 1);
    maxSamps = gazeArray(idx_longest).NumSamples;
    
% tabulate

    % [lx, ly, rx, ry, t, missingLeft, missingRight, absent]
    lx = nan(maxSamps, numGaze);
    ly = nan(maxSamps, numGaze);
    rx = nan(maxSamps, numGaze);
    ry = nan(maxSamps, numGaze);
    t = nan(maxSamps, numGaze);
    missingLeft = false(maxSamps, numGaze);
    missingRight = false(maxSamps, numGaze);
    absent = false(maxSamps, numGaze);
    
    g = nan(maxSamps, 8, numGaze);
    for i = 1:numGaze
    
        s1 = 1;
        s2 = gazeArray(i).NumSamples;
        
        lx(s1:s2, i) = gazeArray(i).LeftX;
        ly(s1:s2, i) = gazeArray(i).LeftY;
        rx(s1:s2, i) = gazeArray(i).RightX;
        ry(s1:s2, i) = gazeArray(i).RightY;
        t(s1:s2, i) = gazeArray(i).Time;
        missingLeft(s1:s2, i) = gazeArray(i).LeftMissing;
        missingRight(s1:s2, i) = gazeArray(i).RightMissing;
        absent(s1:s2, i) = false;
        absent(s1:s2, i) = true;
        
    end
    
    % take the median timestamp value for each sample
    t = nanmedian(t, 2);

% make gaze object

    gaze = etGazeDataBino;
    gaze.Import(lx, ly, rx, ry, t, missingLeft, missingRight, absent);
    
end