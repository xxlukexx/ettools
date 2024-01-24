function imgd = etDilateMaskAOI_image(img, def, radius)

    etAssertAOIDef(def)
    tol = 10;
    
    % split input image into RGB channels
    imgd_r = img(:, :, 1);
    imgd_g = img(:, :, 1);
    imgd_b = img(:, :, 1);
    
    % get width and height
    w = size(img, 2);
    h = size(img, 1);
    
    % loop through each AOI
    numAOIs = size(def, 1);
    for a = 1:numAOIs
        
    % each AOI can be defined by one or more colour values. Process
    % each colour value to produce a binary mask for this AOI
    
        numCol = length(def{a, 2}); 
        
        % blank mask
        mask = false(size(img, 1), size(img, 2));
        
        % loop through colours in def for this AOI and build the mask on
        % each interation
        for c = 1:numCol
            
            % extract colour value from definition
            col = def{a, 2}{c};

            % find AOI pixels for this def
            mask = mask | roiRGB(img, col, tol);
            
        end        
        
        % dilate mask
        se = strel('disk', radius, 0);
        maskd = imdilate(mask, se);
        
        % get x, y coords of all white pixels in the dilated mask
        [y, x] = ind2sub([h, w], find(maskd));

        % colour output image using the mask and AOI def
        imgd_r(maskd) = col(1);
        imgd_g(maskd) = col(2);
        imgd_b(maskd) = col(3);
        
%         subplot(1, 2, 1)
%         imshow(mask)
%         subplot(1, 2, 2)
%         imshow(cat(3, imgd_r, imgd_g, imgd_b))
        
    end

    % combine RGB channels for output image
    imgd = cat(3, imgd_r, imgd_g, imgd_b);


end