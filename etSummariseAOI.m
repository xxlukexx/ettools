function [smry, tab_tall, tab_wide, looks] = etSummariseAOI(varargin)

    % convert in matrix to logical 
    varargin{1} = logical(varargin{1});

    if nargout == 4
        % include looks
        [smry, tab_tall, tab_wide, looks] = etSummariseIn('aoi', varargin{:});
    else
        [smry, tab_tall, tab_wide] = etSummariseIn('aoi', varargin{:});
    end
    
end