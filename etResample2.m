% [mbr, tbr, tr] = ETRESAMPLE(mb, tb, fs_new)
%
% Resamples eye tracking data in the Tobii analytics format (n x 26
% matrix). MB is a main buffer, TB is a time buffer, FS_NEW is the new
% sampling rate. 
%
% Output arguments are MBR, TBR and TR, the resamples main buffer, time
% buffer and time vector, respectively. 
%
% When upsampling, no interpolation is performed, so samples will be
% duplicated. 
%
% When downsampling, missing data is relaced with NaNs, removed, and then
% the mean is taken. If all data are missing for a particular output sample
% then the output sample will also be missing. Otherwise, the output sample
% will be the mean of the available (non-missing) gaze data. 
%
% For each output sample, the minimum gaze validity code is taken from all
% input samples. This means that if any samples were non-missing (validity
% code <4), the validity code of the output sample will be <4. We cannot
% completely recreate the validity codes as the Tobii eye tracker produces
% them, because we cannot access the confidence as to which eye is which.
% For most analyses, it's fine to treat the presence/absence of NaNs in the
% gaze data for each eye as a binary validity code. 
function [mbr, tbr, tr] = etResample2(mb, tb, fs_new)

    mbr = mb;
    tbr = tb;
    tr = [];
    
    if isempty(mb) || isempty(tb)
        return
    end

    % replace missing with NaNs 
    mb = etPreprocess(mb, 'removemissing', true);
    
    % make time vector, total duration
    t = etTimeBuffer2Secs(tb);
    dur = t(end);
    
    % make time vector for resampled data
    tr(:, 1) = 0:1 / fs_new:dur;
    numSamps_new = length(tr);

    % for each element of the resampled data, find the left-hand edge of
    % the corresponding period in the original data
    s1 = arrayfun(@(x) find(t >= x, 1, 'first'), tr);
    s2 = [s1(2:end) - 1; length(t)];
    s2(s2 < s1) = s1(s2 < s1);

    % for all columns except sample validity, take the mean value inside
    % each period of original data. 
    mbr = nan(numSamps_new, 26);
    for s = 1:numSamps_new
        
        % left eye 
        mbr(s, 1:12) = nanmean(mb(s1(s):s2(s), 1:12), 1);
        mbr(s, 13) =  min(mb(s1(s):s2(s), 13), [], 1);
        
        % right eye
        mbr(s, 14:25) = nanmean(mb(s1(s):s2(s), 14:25), 1);
        mbr(s, 26) =  min(mb(s1(s):s2(s), 26), [], 1);
        
    end
    
    % recreate time buffer using output sample indices
    tbr = tb(s1, :);

end