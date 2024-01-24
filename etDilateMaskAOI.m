function etDilateMaskAOI(path_aoi, def, radius)

% check input args

    % check def 
    etAssertAOIDef(def)

    % check AOI path
    if ~exist(path_aoi, 'file')
        error('File not found: %s', path_aoi)
    end

% check for valid formats, and determine if video of image

    valFormats_img = {'png', 'tiff', 'jpeg', 'jpg'};
    valFormats_vid = {'avi', 'mp4', 'mpeg4', 'mov'};
    
    % get extension
    [pth, fil, ext] = fileparts(path_aoi);
    
    % compare extension (stripping off leading '.') to valid formats
    isImg = ismember(ext(2:end), valFormats_img);
    isVid = ismember(ext(2:end), valFormats_vid);
    
    if ~isImg && ~isVid
        error('Unsupported file format.')
    end
    
% make output filename

    % append '_dilated' to filename
    fil_out = sprintf('%s_dilated%s', fil, ext);
    if isVid, ext = 'mp4'; end
    path_out = fullfile(pth, fil_out);
    
% call relevant function depending upon format

    if isImg
        
        % attempt to load
        try
            img = imread(path_aoi);
        catch ERR
            error('Error reading image file: %s', ERR.message)
        end
        
        imgd = etDilateMaskAOI_image(img, def, radius);
        imwrite(imgd, path_out)
        
        % report result and plot original and dilated AOI
        fprintf('Wrote dilated image file to: %s\n', path_out);
        figure
        subplot(1, 2, 1)
        imshow(img)
        title('Original')
        subplot(1, 2, 2)
        imshow(imgd)
        title('Dilated')
        
    elseif isVid
        
        % attempt to load
        try
            vr = VideoReader(path_aoi);
        catch ERR
            error('Error reading video file: %s', ERR.message)
        end
        
        % make output video
        try
            vw = VideoWriter(path_out, 'MPEG-4');
            vw.FrameRate = vr.FrameRate;
            vw.Quality = 90;
            open(vw)
        catch ERR
            error('Error writing video file: %s', ERR.message)
        end
        
        % loop through and dilate each frame
        fr = 1;
        while vr.hasFrame
            
            frame = vr.readFrame;
            frame_dl = etDilateMaskAOI_image(frame, def, radius);
            vw.writeVideo(frame_dl);
            fr = fr + 1;
            
            if mod(fr, 20) == 0
                fprintf('Frame %d...\n', fr);
            end
            
        end
        
        close(vw)
            
    end

end