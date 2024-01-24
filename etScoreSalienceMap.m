function in = etScoreSalienceMap(gaze, img)

% general setup
    
    % check gaze is a teGazeData instance
    if ~exist('gaze', 'var') || ~isa(gaze, 'etGazeData') || isempty(gaze)
        error('Must supply ''gaze'', a non-empty teGazeData instance.')
    end
    
    if isempty(gaze)
        in = [];
        return
    end
    
    % if the passed map is RGB, convert to grey
    if size(img, 3) ~= 1
        img = rgb2gray(img);
    end
    
    % currently only supports one map - this may change in future (to allow
    % multiple salience features to be nested within one object - e.g.
    % contrast, colour)
    numMaps = 1;
    in = nan(gaze.NumSamples, gaze.NumSubjects, numMaps);

    % pull time, coords and image data into local vars
    
        % get [x, y] coords
        x = gaze.X;
        y = gaze.Y;

        % get image data
        w = size(img, 2);
        h = size(img, 1);

    % get x, and y values

        x = 1 + floor(x .* w);
        y = 1 + floor(y .* h);

        % put x any y into tall vectors, to speed up pixel lookup by
        % avoiding a loop
        x = reshape(x, [], 1);
        y = reshape(y, [], 1);
        missing = reshape(gaze.Missing | gaze.Absent, [], 1);

    % look up intensity (salience) values at each gaze point

        % convert [x, y] coords to single indices within the AOI image
        % data
        idx = sub2ind([h, w], y, x);

        % find nans (missing gaze samples) and replace with ones. We
        % will use the index of nan values later on to remove these
        % zeros and put nans back in their place
        idx(missing) = 1;

        % look up the pixel value at each index (gaze point)
        intensity = double(img(idx));

        % put missing data samples back in as nans
        intensity(gaze.Missing) = false;

        % convert back from single column vector to samples x
        % sub matrix
        in = reshape(intensity, gaze.NumSamples, gaze.NumSubjects);

end