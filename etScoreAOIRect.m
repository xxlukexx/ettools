function in = etScoreAOIRect(gaze, rect)
    
% general setup

    % check format of rect argument - must be a [n x 4] vector/matrix,
    % where n is the number of AOIs, and columns represent [x1, x2, y1, y2]
    if ~isnumeric(rect) || size(rect, 2) ~= 4
        error('''rect'' argument must be an [n x 4] vector or matrix, with a row for each AOI.')
    end
    
    % check gaze is a teGazeData instance
    if ~exist('gaze', 'var') || ~isa(gaze, 'etGazeData') || isempty(gaze)
        error('Must supply ''gaze'', a non-empty teGazeData instance.')
    end
    
    if isempty(gaze)
        in = [];
        return
    end
    
    % determine number of AOIs
    numAOIs = size(rect, 1);

% score AOI(s)

    % find gaze samples inside the AOI. Samples from the left, right or the
    % averaged gaze position are counted if they are in the AOI
    in = nan(gaze.NumSamples, gaze.NumSubjects, numAOIs);
    
    for a = 1:numAOIs
        
        % get gaze samples within AOI
        in(:, :, a) = (...    
                gaze.X >= rect(a, 1) &...                                            % average gaze 
                gaze.X <= rect(a, 3) &...
                gaze.Y >= rect(a, 2) &...
                gaze.Y <= rect(a, 4)...
            ) |...
            (...
                gaze.LeftX >= rect(a, 1) &...                                           % left eye
                gaze.LeftX <= rect(a, 3) &...
                gaze.LeftY >= rect(a, 2) &...
                gaze.LeftY <= rect(a, 4)...
            ) |...
            (...
                gaze.RightX >= rect(a, 1) &...                                           % right eye
                gaze.RightX <= rect(a, 3) &...
                gaze.RightY >= rect(a, 2) &...
                gaze.RightY <= rect(a, 4)...
            );
    end
    
    in = logical(in);
    
end
    
