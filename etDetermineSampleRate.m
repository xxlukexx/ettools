function [sampleRate, secondsPerSample] = etDetermineSampleRate(timeBuffer)
%ETDETERMINESAMPLERATE Estimate rate while excluding gaps/outlier deltas.
if isempty(timeBuffer)
    warning('ET:EmptyTimeBuffer', 'Time buffer is empty.')
    sampleRate = nan;
    secondsPerSample = nan;
    return
end
timestamps = double(timeBuffer(:, 1));
delta = [nan; diff(timestamps)];
jumpCriterionUs = 100 * 1000;
outlier = detectOutliers(delta, 5) | delta > jumpCriterionUs;
delta(outlier) = [];
secondsPerSample = nanmean(delta) / 1e6;
sampleRate = 1 / secondsPerSample;
end
