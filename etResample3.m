function [mbr, tbr, tr] = etResample3(mb, tb, fs_new)
%ETRESAMPLE3 Resample a Tobii Analytics-style 26-column gaze buffer.
%
% [mbr, tbr, tr] = etResample3(mb, tb, fs_new)
%
% Continuous columns are averaged within regular output bins after missing
% gaze samples have been converted to NaNs. Per-eye validity columns 13 and
% 26 use the minimum non-NaN value in each bin. Empty bins use the nearest
% original sample, which also provides the corresponding output time-buffer
% row. This supports both downsampling and non-interpolating upsampling.

    mbr = mb;
    tbr = tb;
    tr = [];

    if nargin < 3
        error('etResample3 requires mb, tb, fs_new.');
    end

    if isempty(mb) || isempty(tb)
        mbr = [];
        tbr = [];
        return
    end

    if ~isnumeric(mb) || size(mb, 2) < 26
        error('mb must be a numeric matrix with at least 26 columns.');
    end

    if ~isscalar(fs_new) || ~isfinite(fs_new) || fs_new <= 0
        error('fs_new must be a positive finite scalar.');
    end

    % Ignore acquisition-specific columns beyond the Tobii Analytics schema.
    mb = mb(:, 1:26);
    mb = etPreprocess(mb, 'removemissing', true);

    t = double(etTimeBuffer2Secs(tb));
    t = t(:);

    if isempty(t) || numel(t) ~= size(mb, 1)
        error('Time vector length must match number of rows in mb.');
    end

    if any(~isfinite(t))
        error('Time vector contains non-finite values.');
    end

    if any(diff(t) < 0)
        error('Time vector is not monotonic nondecreasing.');
    end

    t0 = t(1);
    t1 = t(end);
    dt = 1 / fs_new;

    if t1 <= t0
        tr = t0;
        mbr = mb(1, :);
        tbr = tb(1, :);
        return
    end

    numSamps_new = floor((t1 - t0) * fs_new) + 1;
    tr = t0 + (0:numSamps_new - 1)' * dt;
    tr = tr(tr <= t1 + dt / 2);
    numSamps_new = numel(tr);

    mbr = nan(numSamps_new, 26);
    idxNearest = nan(numSamps_new, 1);
    startIdx = 1;
    nIn = numel(t);

    for s = 1:numSamps_new
        binStart = tr(s);

        if s < numSamps_new
            binEnd = tr(s + 1);
            idx = find(t(startIdx:end) >= binStart & ...
                t(startIdx:end) < binEnd) + startIdx - 1;
        else
            idx = find(t(startIdx:end) >= binStart) + startIdx - 1;
        end

        if isempty(idx)
            [~, idx] = min(abs(t - binStart));
        else
            startIdx = idx(1);
        end

        idxNearest(s) = idx(1);
        thisMb = mb(idx, :);

        mbr(s, 1:12) = localNanMean(thisMb(:, 1:12));
        mbr(s, 13) = localNanMin(thisMb(:, 13));
        mbr(s, 14:25) = localNanMean(thisMb(:, 14:25));
        mbr(s, 26) = localNanMin(thisMb(:, 26));
    end

    if size(tb, 1) ~= nIn
        error('tb must have same number of rows as mb.');
    end

    idxNearest = max(1, min(nIn, idxNearest));
    tbr = tb(idxNearest, :);
end

function out = localNanMean(x)

    out = nan(1, size(x, 2));
    for c = 1:size(x, 2)
        values = x(:, c);
        values = values(~isnan(values));
        if ~isempty(values)
            out(c) = mean(values);
        end
    end
end

function out = localNanMin(x)

    x = x(:);
    x = x(~isnan(x));
    if isempty(x)
        out = nan;
    else
        out = min(x);
    end
end
