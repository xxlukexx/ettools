function etAssertAOIDef(val)

    % check correctness. The definition should be a [n x 2] cell
    % array, where n is the number of AOIs. The first column is a
    % cell array of strings listing the AOI names, and the second
    % column is a cell array of cell arrays, each containing one or
    % more [1 x 3] RGB colour values, representing the colour(s) in
    % that AOI
    if ~iscell(val) || ~size(val, 2) == 2 ||...
            ~iscellstr(val(:, 1)) ||...
            (~all(cellfun(@isnumeric, val(:, 2))) &&...
            ~all(cellfun(@iscell, val(:, 2))))
        error('AOIDefinition must be a [n x 2] cell array. The first column is a list of aoi names, and the second a list of colours.')
    end
            
end