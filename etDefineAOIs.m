function def = etDefineAOIs(img, def, startTime)

    if ~exist('def', 'var') || isempty(def)
        def = {};
    end
    
    if ~exist('startTime', 'var') || isempty(startTime)
        startTime = 0;
    end
    
    if ischar(img)
        if ~exist(img, 'file'), error('File not found.'), end
        % try to determine type by extension
        [~, ~, ext] = fileparts(img);
        switch lower(ext(2:end))
            case {'png', 'jpg', 'jpeg', 'bmp', 'gif'}
                try
                    img = imread(img);
                catch ERR
                    error('Could not load image.')
                end
            case {'mp4', 'mov', 'avi'}
                try
                    vr = VideoReader(img);
                    vr.CurrentTime = startTime;
                    img = readFrame(vr);
                catch ERR
                    error('Could not load movie.')
                end
            otherwise
                error('File formats supported are images (png, jpg, bmp, gif) or movies (mp4, mov, avi).')
        end
    end
               
    % get AOI name and feature nampwd
    aoiFeatureName = input('Enter AOI/feature name (e.g. "aoi_face"): >', 's');

    % get number of colours
    numColHappy = false;
    while ~numColHappy
        aoiNumCols =...
            str2double(input('How many colours in this AOI?', 's'));
        numColHappy = ~isempty(aoiNumCols);
    end
    
    fig = figure('Name', 'Pick Colour', 'menubar', 'none');
    imshow(img);
    
    [x, y] = ginput(aoiNumCols);
    
    defCol = cell(1, aoiNumCols);
    mask = false(size(img, 1), size(img, 2));
    for c = 1:aoiNumCols
        defCol{c} = shiftdim(img(round(y(c)), round(x(c)), :), 1);
        mask = mask | roiRGB(img, defCol{c});
    end
    
    def{end + 1, 1} = aoiFeatureName;
    def{end, 2} = defCol;
    
    test = img;
    mask = repmat(mask, 1, 1, 3);
    test(~mask) = 0;
    imshow(test);
    
end