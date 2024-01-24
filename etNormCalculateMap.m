function [map, scores] = etNormCalculateMap(gaze, resolution)

    % default resolution is 20 x 20 
    if ~exist('resolution', 'var') || isempty(resolution)
        res_x = linspace(0, 1, 21);
        res_y = linspace(0, 1, 21);
    else
        if ~isvector(resolution) && length(resolution) ~= 2
            error('resolution must be a two-element [w, h] vector.')
        else
            res_x = resolution(1);
            res_y = resolution(2);
        end
    end
    
% create gaze probability map
      
    % preallocate map storage
    map = nan(length(res_x) - 1, length(res_y) - 1, gaze.NumSamples);
    
    % calculate sample-by-sample maps
    for s = 1:gaze.NumSamples
        x = gaze.X(s, :)';
        y = gaze.Y(s, :)';
        map(:, :, s) =...
            histcounts2(y, x, res_y, res_x, 'Normalization', 'probability');
    end
    
% score each individual against the map

    if nargout == 2
        scores = etNormScoreSubjectsAgainstMap(gaze, map);
    end

end