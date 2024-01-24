function etPlotAOIEffectSizes(compare, measure, img_aoi, tall, wide, def, img_stim)

%     measure = 'peakLook';
    
% setup and load

    % default to prop looking
    if ~exist('measure', 'var') || isempty(measure)
        measure = 'propInAOI';
    end

    % only plot stimulus if a path has been supplied
    plotStim = exist('img_stim', 'var') && ~isempty(img_stim);
    if plotStim
        img_stim = checkloadimage(img_stim);
        % convert to grayscale
        img_stim =  rgb2gray(img_stim);
    end
    
    % load AOI image
    img_aoi = checkloadimage(img_aoi);
    
    % resize stim if needed
    if ~isequal(size(img_aoi), size(img_stim))
        img_stim = imresize(img_stim, [size(img_aoi, 1), size(img_aoi, 2)]);
    end
    
    % generate AOI alpha
    r = img_aoi(:, :, 1);
    g = img_aoi(:, :, 2);
    b = img_aoi(:, :, 3);
    alpha_aoi = r == 0 & g == 0 & b == 0;
    
% t-tests and ES

    % check there are two levels of compare
    if length(unique(tall.(compare))) ~= 2
        warning('Cannot perform t-test with more than two levels of the compare (%s) variable.',...
            compare)
    else

        numAOIs = size(def, 1);
        res = struct;
        for a = 1:numAOIs
            
            res(a).aoi = def(a, 1);

            % filter for this AOI
            idx_aoi = strcmpi(tall.aoi, def{a, 1});
            tmp = tall(idx_aoi, :);

            % get subscripts for levels of compare variable
            [comp_u, ~, comp_s] = unique(tmp.(compare));

            % t-test
            [~, p, ci, stats] = ttest2(tmp.(measure)(comp_s == 1),...
                tmp.(measure)(comp_s == 2));
            res(a).p = p;
            res(a).ci = ci;
            res(a).t = stats.tstat;
            res(a).df = stats.df;
            res(a).sd = stats.sd;
            
            % cohen's D
            mu = accumarray(comp_s, tmp.(measure), [], @nanmean);
            x1 = tmp.(measure)(comp_s == 1);
            x2 = tmp.(measure)(comp_s == 2);
            [~, sdp] = sdpooled(x1, x2);
            res(a).meandiff = diff(mu);
            res(a).d = abs(diff(mu) / sdp);        
        
        end
        
    end
    
% find AOI centroids

    cx = nan(numAOIs, 1);
    cy = nan(numAOIs, 1);
    for a = 1:numAOIs

        % mask just this AOI 
        numAOICols = length(def(a, 2));
        mask = false(size(img_aoi, 1), size(img_aoi, 2));
        for c = 1:numAOICols
            mask = mask | roiRGB(img_aoi, def{a, 2}{c});
        end
        
        % find centroid
        [cx(a), cy(a)] = etFindAOICentroid(mask);

    end
    
    % rescale centroid from normalised to pixels
    w = size(img_aoi, 2);
    h = size(img_aoi, 1);
    cx = cx .* w;
    cy = cy .* h;
        
    
% plot

    img_bar = LEAP_ET_ns_plotBar(wide, measure, 'diag', [], [], [], def(:, 1));
    fig = figure('WindowStyle', 'docked');
    
    subplot(1, 2, 1)
    
        % optionally plot stim
        if plotStim, imshow(img_stim), end

        % plot AOI
        hold on
        h = imshow(img_aoi, 'border', 'tight');
        h.AlphaData = (alpha_aoi * .7) + (~alpha_aoi * .4);
        
        % plot results over each AOI
        for a = 1:numAOIs
            
            if res(a).meandiff > 0
                strGrp = sprintf('%s>%s', comp_u{2}, comp_u{1});
            elseif res(a).meandiff < 0 
                strGrp = sprintf('%s>%s', comp_u{1}, comp_u{2});
%             else
%                 strGrp = sprintf('%s~%s', comp_u{1}, comp_u{2});
            end
            str = sprintf('%s\n%s\nd=%.2f [p=%.3f]', strGrp, def{a, 1}, res(a).d, res(a).p);
            h1 = text(cx(a) + 1, cy(a) + 1, str, 'Color', 'k', 'FontSize', 20, 'HorizontalAlignment', 'center');
            h2 = text(cx(a), cy(a), str, 'Color', 'w', 'FontSize', 20, 'HorizontalAlignment', 'center');
            if res(a).p < .05
                h1.FontWeight = 'bold';
                h2.FontWeight = 'bold';
            end
%             pos = textBounds(str, gca);
%             rectangle('Position', pos, 'FaceColor', [0, 0, 0, .7], 'LineWidth', 2, 'EdgeColor', 'w')
            
        end
        
    subplot(1, 2, 2)
    
        imshow(img_bar, 'border', 'tight')

end