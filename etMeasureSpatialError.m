function [acc, prec] = etMeasureSpatialError(gaze, kx, ky, varargin)
% measures accuracy and precision in an etGazeData instance. If gaze is a cell
% array of etGazeData instances then the function will call itself on each
% element. 
%
% [kx, ky] is the spatial location in gaze coords of a known stimulus.
% Accuracy and precision will be measured relative to this point.
%
% acc is the RMS 2D euclidean distance from [kx, ky] to each gaze point. 
%
% prec is the RMS 2D euclidean distance from the centroid of all gaze
% points. 
%
% Missing data is removed. If this removes all samples, NaNs are returned
% for acc and prec.
%
% Optional input args:
%
%       t1 - allows measuring spatial error on a subset of the data,
%       defined as samples between t1 and t2. If t1 is specified, t2 is
%       required. 
%
%       t2 - as above
%
%       units - causes error to be calculated using units other than the
%       default normalised gaze coords ('norm'). Options are 'deg' (degrees of
%       visual angle) or 'screen' (screen units). Note that for either
%       'deg', or 'screen' to return non-NaN, the etGazeData's
%       .ScreenDimensions vector must be set. Further, for 'deg' to work,
%       the .DistanceFromScreen property must be set. Note that [kx, ky]
%       must be in the same units.
    
    parser      =   inputParser;
    checkTime   =   @(x) isnumeric(x) && isscalar(x) && x >= 0;
    checkUnits  =   @(x) ismember(x, {'norm', 'deg', 'screen'});
    addParameter(   parser, 't1',                   [],         checkTime       )
    addParameter(   parser, 't2',                   [],         checkTime       )    
    addParameter(   parser, 'units',                'norm',     checkUnits      )
    parse(          parser, varargin{:});
    t1          =   parser.Results.t1;
    t2          =   parser.Results.t2;
    units       =   parser.Results.units;

    % if gaze is cell array, call this function for each element
    if iscell(gaze)
        [acc, prec] =...
            cellfun(@(x) etMeasureSpatialError(x, varargin{:}), gaze);
        return
    end

    % check gaze format
    if ~isa(gaze, 'etGazeData')
        error('gaze must be a teGazeData instance (mono or bino).')
    end

    % check that if either t1/t2 are specified, the other is also
    if (~isempty(t1) && isempty(t2)) || (isempty(t1) && ~isempty(t2))
        error('If t1 if specified, t2 must be (and vice versa).')
    end
    
    % if t1/t2 not specified, use the entire teGazeData length
    if isempty(t1)
        t1 = 0;
    end
    if isempty(t2)
        t2 = gaze.Time(end);
    end    
    
    % check t1/t2 bounds
    if ~isempty(t1)
        if t1 > gaze.Time(end) 
            warning('t1 (%.3f) is greater than gaze data length (%.3f).',...
                t1, gaze.Time(end))
        elseif t2 > gaze.Time(end)
            warning('t2 (%.3f) is greater than gaze data length (%.3f). Will use gaze data length instead of t2.',...
                t2, gaze.Time(end))
            t2 = gaze.Time(end);
        end
    end
    
    % get x, y
    idx = gaze.Time >= t1 & gaze.Time <= t2;
    switch units
        case 'norm'
            x = gaze.X(idx);
            y = gaze.Y(idx);
        case 'screen'
            x = gaze.XScreen(idx);
            y = gaze.YScreen(idx);
        case 'deg'
            x = gaze.XDeg(idx);
            y = gaze.YDeg(idx);
    end
    missing = gaze.Missing(idx) | gaze.Absent(idx);
    if all(missing)
        acc = nan;
        prec = nan;
        return
    end
    
    % remove missing
    x(missing) = [];
    y(missing) = [];
    
    % for each sample, calculate euclidean distance from known x, y
    dxk = x - kx;
    dyk = y - ky;
    dis_k = sqrt((dxk .^ 2) + (dyk .^ 2));
    
    % accuracy
    acc = rms(dis_k);
    
    % precision
    cx = mean(x);
    cy = mean(y);
    dxg = x - cx;
    dyg = y - cy;
    dis_g = sqrt((dxg .^ 2) + (dyg .^ 2));
    prec = rms(dis_g);

end