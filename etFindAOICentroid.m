function [cx, cy] = etFindAOICentroid(bin)

    if ~islogical(bin)
        error('AOI must be binarised - use etBinariseAOI first on RGB images.')
    end
    
    numAOIs = size(bin, 3);
    w = size(bin, 2);
    h = size(bin, 1);
    
    cx = nan(numAOIs, 1);
    cy = nan(numAOIs, 1);
    for a = 1:numAOIs
        
        [y, x] = ndgrid(1:h, 1:w);
        centroid =...
            round(mean([x(logical(bin(:, :, a))), y(logical(bin(:, :, a)))]));        
        cx(a) = centroid(1) / w;
        cy(a) = centroid(2) / h;
        
    end

end