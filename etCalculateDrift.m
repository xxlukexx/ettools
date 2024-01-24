function [accuracy, precision] = etCalculateDrift(pt_x, pt_y, gaze_x, gaze_y)

%     % check input args
%     if ~isnumeric(pt_x) || ~isnumeric(pt_y) || ~isscalar(pt_x) ||...
%             ~isscalar(pt_y)
%         error('pt_x and pt_y must be numeric scalars.')
%     end
%     if ~isnumeric(gaze_x) || ~isnumeric(gaze_y) || ~isvector(gaze_x) ||...
%             ~isvector(gaze_y)
%         error('gaze_x and gaze_y must be numeric vectors.')
%     end

    % calculate euclidean distance from each gaze point to the calibration
    % point - this is the accuracy of each sample
    distanceToPoint = sqrt(((pt_x - gaze_x) .^ 2) + ((pt_y - gaze_y) .^ 2));
    
    % for each calibration point separately, calculate the centroid of all
    % gaze points around it
    sig = arrayfun(@(x, y) sprintf('%.30f_%.30f', x, y), pt_x, pt_y,...
        'UniformOutput', false);
    [pt_u, ~, pt_s] = unique(sig);
    numPoints = length(pt_u);
    numSamps = length(pt_x);
    cx = nan(numSamps, 1);
    cy = nan(numSamps, 1);
    for p = 1:numPoints
        idx = pt_s == p;
        cx(idx) = mean(gaze_x(idx));
        cy(idx) = mean(gaze_y(idx));
    end
    
    % caluclate distance from each gaze sample to the centroid of all gaze
    % samples
%     cx = mean(gaze_x);
%     cy = mean(gaze_y);
    distanceToCentroid = sqrt(((gaze_x - cx) .^ 2) + ((gaze_y - cy) .^ 2));
    
    % calculate accuracy (mean spatial error of all gaze points)
    accuracy = nanmean(distanceToPoint);
    
    % calculate precision (rms of all gaze points)
    precision = nanrms(distanceToCentroid);
    
end