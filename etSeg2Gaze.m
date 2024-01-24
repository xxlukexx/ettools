function gaze = etSeg2Gaze(seg, dur)

    if ~exist('dur', 'var') 
        dur = [];
    end

    [t, lx, ly, rx, ry, missingLeft, missingRight, absent, id, wave] =...
        etTabulateSeg(seg, dur);
%     sig_id = cellfun(@(id, wave) sprintf('%s_%s', id, wave), id, wave,...
%         'UniformOutput', false);
    gaze = etGazeDataBino;
    gaze.Import(lx, ly, rx, ry, t, missingLeft, missingRight, absent);

end 