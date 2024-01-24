function dist = etScoreAOIDistanceMask(in, gaze, img, def, postInterpX,...
    postInterpY)
% find centroid of AOI, and calculate distance between each gaze point
% and centroid for a measure of distance

    % optionally use post-interpolation [x, y] coords. This allows for
    % interpolated gaze coords to be used so that distance estimates are in
    % line with AOI scores. If not supplied, simply use the raw [x, y] from
    % the etGazeData instance
    if ~exist('postInterpX', 'var') || isempty(postInterpX)
        x = gaze.X;
    else 
        x = postInterpX;
    end
    
    if ~exist('postInterpY', 'var') || isempty(postInterpY)
        y = gaze.Y;
    else
        y = postInterpY;
    end

    % convert normalised gaze coords to pixels based on image size
    w = size(img, 2);
    h = size(img, 1);
    x = reshape(x, [], 1);
    y = reshape(y, [], 1);
        
    % binarise image
    img_bin = etBinariseAOIMask(img, def);

    % find centroid, convert to normalised coords
    [cx, cy] = etFindAOICentroid(img_bin);

    % calculate distance from each gaze point, to each AOI centroid
    numAOIs = size(def, 1);
    dist = nan(gaze.NumSamples, gaze.NumSubjects, numAOIs);
    for a = 1:numAOIs
        
        idx = reshape(in(:, :, a), [], 1);
        dx = x;
        dy = y;
        dx(~idx) = nan;
        dy(~idx) = nan;
        dx = dx - cx(a);
        dy = dy - cy(a);
        dist_tmp = sqrt((dx .^ 2) + (dy .^ 2));
        dist(:, :, a) = reshape(dist_tmp, gaze.NumSamples, []);
    end

end