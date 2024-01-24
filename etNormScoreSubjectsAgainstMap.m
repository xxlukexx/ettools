function scores = etNormScoreSubjectsAgainstMap(gaze, map)

    % calculate resolution of map
    res_x = size(map, 2);
    res_y = size(map, 1);

    scores = nan(gaze.NumSamples, gaze.NumSubjects);
    for s = 1:gaze.NumSamples
        
        % round gaze to nearest map bin
        x = floor(gaze.X(s, :) * res_x) + 1;
        y = floor(gaze.Y(s, :) * res_y) + 1;  
        
        % put x any y into tall vectors, to speed up pixel lookup by
        % avoiding a loop
        x = reshape(x, [], 1);
        y = reshape(y, [], 1);
        missing = reshape(gaze.Missing(s, :) | gaze.Absent(s, :), [], 1);

    % look up colour values at each gaze point

        % convert [x, y] coords to single indices within the AOI image
        % data
        idx = sub2ind([res_y, res_x], y, x);

        % find nans (missing gaze samples) and replace with ones. We
        % will use the index of nan values later on to remove these
        % zeros and put nans back in their place
        idx(missing) = 1;
        
        % filter map for just this sample
        map_tmp = map(:, :, s);
        
        % look up probability scores for each subject on the map
        scores(s, :) = map_tmp(idx);
        
        % replace NaNs
        scores(s, missing) = nan;

        
    end



end