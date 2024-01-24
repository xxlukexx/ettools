function etWriteDataDrivenAOIsToDisk(aoi, def, path_out, label)

    if ~exist('label', 'var') || isempty(label)
        label = sprintf('DDAOI_%s', datestr(now, 30));
    end
    
    tryToMakePath(path_out);
    filename_aoi = sprintf('daoi_%s.png', label);
    filename_def = sprintf('aoidef_%s.mat', label);
    file_aoi = fullfile(path_out, filename_aoi);
    file_def = fullfile(path_out, filename_def);
    
    imwrite(aoi, file_aoi)
    save(file_def, 'def')
    
    fprintf('Wrote DD-AOI image in PNG format to: %s\n', file_aoi);
    fprintf('Wrote DD-AOI definition in MAT format to: %s\n', file_def);
    
end