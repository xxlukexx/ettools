function bin = etBinariseAOIMask(mask, def)

%%

% tic

    % if the mask has an alpha channel, remove it (we only work with RGB)
    if size(mask, 3) == 4
        mask = mask(:, :, 1:3);
    end

    numAOIs = size(def, 1);
    w = size(mask, 2);
    h = size(mask, 1);
    bin = false(h, w, numAOIs);
    tol = 10;
    
    mask = double(mask);
%     mask_r = mask(:, :, 1);
%     mask_g = mask(:, :, 2);
%     mask_b = mask(:, :, 3);
        
    for a = 1:numAOIs

        % determine number of colours in this AOI
        numCols = size(def{a, 2}, 2);
        idx_aoi = false(h, w, numCols);
        for c = 1:numCols

            % pull RGB values from the def
            def_r = double(def{a, 2}{c}(1));
            def_g = double(def{a, 2}{c}(2));
            def_b = double(def{a, 2}{c}(3));
            def_col = cat(3, def_r, def_g, def_b);
            
%             % compare against AOI pixel values
%             idx = ...
%                 abs(mask_r - def_r) < tol &...
%                 abs(mask_g - def_g) < tol &...
%                 abs(mask_b - def_b) < tol;
                        
%             idx = all(abs(mask - def_col) < tol, 3);
            idx_aoi(:, :, c) = all(abs(mask - def_col) < tol, 3);
                        
%             bin(:, :, a) = bin(:, :, a) | idx;
            
        end
        
        bin(:, :, a) = any(idx_aoi, 3);

    end    
    
% toc

%%

end