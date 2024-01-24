function [in, cols, dist] = etScoreAOIMask(gaze, img, def, colourTolerance)
    
    % check gaze is a teGazeData instance
    if ~exist('gaze', 'var') || ~isa(gaze, 'etGazeData') || isempty(gaze)
        error('Must supply ''gaze'', a non-empty teGazeData instance.')
    end
    
    if isempty(gaze)
        in = [];
        return
    end
    
    % colourTolerance is the tolerance in pixel values (out of 255) that
    % will count as a 'hit' in the AOI. This is used when an image format
    % doesn't represent colours accurately (e.g. mp4) or where the colour
    % space has changed. Default to 1 - only score pixels that exactly
    % match the colour in the AOI definition
    if ~exist('colourTolerance', 'var') || isempty(colourTolerance)
        colourTolerance = 1;
    end
    
    % determine number of AOIs
    numAOIs = size(def, 1);
    in = nan(gaze.NumSamples, gaze.NumSubjects, numAOIs);

    % pull time, coords and image data into local vars
    
        % get [x, y] coords
        gx = gaze.X;
        gy = gaze.Y;

        % get image data
        w = size(img, 2);
        h = size(img, 1);

        % split image into separate variables for each colour channel.
        % This allows us to use one set of gaze point indices to look
        % up the intensity value in each channel
        img_r = img(:, :, 1);
        img_g = img(:, :, 2);
        img_b = img(:, :, 3);

    % get x, and y values

        x = 1 + floor(gx .* w);
        y = 1 + floor(gy .* h);

        % put x any y into tall vectors, to speed up pixel lookup by
        % avoiding a loop
        x = reshape(x, [], 1);
        y = reshape(y, [], 1);
        missing = reshape(gaze.Missing | gaze.Absent, [], 1);

    % look up colour values at each gaze point

        % convert [x, y] coords to single indices within the AOI image
        % data
        idx = sub2ind([h, w], y, x);

        % find nans (missing gaze samples) and replace with ones. We
        % will use the index of nan values later on to remove these
        % zeros and put nans back in their place
        idx(missing) = 1;

        % for each rgb channel, look up the pixel value at each index
        % (gaze point)
        r = double(img_r(idx));
        g = double(img_g(idx));
        b = double(img_b(idx));

    % compare colour values against AOI def. AOIs can be comprised of
    % more than one colour, in which case they will have multiple rgb
    % values. So we loop through each AOI and classify gaze points
    % belonging to that AOI.

        % preallocate storage for colour values. The third dimension
        % represents r, g, b colour channels
        cols = repmat(070, size(idx, 1), 3);

        for a = 1:numAOIs

            % determine number of colours in this AOI
            numCols = size(def{a, 2}, 2);
            for c = 1:numCols

                % pull RGB values from the def
                def_r = double(def{a, 2}{c}(1));
                def_g = double(def{a, 2}{c}(2));
                def_b = double(def{a, 2}{c}(3));

                % compare against AOI pixel values
                in_tall = ...
                    abs(r - def_r) < colourTolerance &...
                    abs(g - def_g) < colourTolerance &...
                    abs(b - def_b) < colourTolerance;

                % assign the colour of this AOI to any samples that are
                % in it
                cols(in_tall, :) =...
                    repmat([def_r, def_g, def_b], sum(in_tall), 1);

                % convert to double, to allow storing nans for mising
                % data (logicals can't have nans in them)
                in_tall = double(in_tall);

                % put missing data samples back in as nans
                in_tall(gaze.Missing) = false;

                % convert back from single column vector to samples x
                % sub matrix
                in_mat = reshape(in_tall, gaze.NumSamples, gaze.NumSubjects);

                % update in variable with values from current colour,
                % for current AOI. If we're on the first (or only)
                % colour, do this straight, otherwise use a logical OR
                % to ensure that gaze in ANY aoi colour counts as in
                % that AOI
                if c == 1
                    in(:, :, a) = in_mat;
                else
                    in(:, :, a) = in(:, :, a) == 1 | in_mat == 1;
                end

            end

        end

    % reshape tall colour values into a [num samps x num subs x 3]
    % matrix
    cols = reshape(cols, gaze.NumSamples, gaze.NumSubjects, 3);
    in = logical(in);    
    
%     anyIn = any(shiftdim(in, 2)', 2);
% %     anyIn = in(:, :, 3);
%     cols_plot = (shiftdim(cols, 2)') ./ 255;
%     figure('visible', 'off')
%     imshow(img)
%     hold on
%     scatter(x, y, 50)
%     scatter(x(anyIn), y(anyIn), 100)
%     file_out = fullfile('/users/luke/desktop/aoiplots', sprintf('%03d_%s.png', s, datestr(now, 30)));
%     fastSaveFig(file_out)
%     delete(gcf)

end