function [mb, tb, info] = etPreprocess(mb, varargin)

%mb = ETPREPROCESS applies various preprocessing steps to eye tracking
%data, in the ECK/Task Engine 'mainbuffer' format. Preprocessing steps are
%applied by optional value pair arguments. If multiple arguments are
%specified, the operations will be applied in the order listed below. 
%
%   removemissing       -   Replace samples of missing data, which have a
%                           default value of -1, with NaN. This operation
%                           removes gaze coords, head position, and pupil
%                           data. 
%
%   removeoffscreen     -   Replace samples of gaze which is off-screen with
%                           NaN. Note this does not operate on missing
%                           data, only valid data whose point of gaze is
%                           off-screen. This operation only removes gaze
%                           positions and pupil data, NOT head position. 
%
%   binocularonly       -   Replaces monocular samples of gaze with NaN,
%                           ensuring that only binocular samples are left.
%
%   averageeyes         -   During binocular tracking, averages the gaze
%                           positions of the left and right eyes. During
%                           monocular tracking, copies the gaze positions
%                           from the available eye to the unavailable eye.
%                           If BINOCULARONLY has been specified,
%                           AVERAGEEYES will only average binocular
%                           samples (monocular samples having already been
%                           removed). 
%
%   resampleFs          -   Resample the gaze, head and pupil data to
%                           resampleFs Hz. Must pass a timebuffer as a
%                           separate argument (see below). See also
%                           removegaps (below). Resampling first median
%                           filters the gaze, head pos and pupil data with
%                           a span equal to one sample at the new sampling
%                           rate (so if 120Hz data is being downsampled to
%                           60Hz, the 120Hz will first be median filtered
%                           with a span of 2 samples). Then a new vector of
%                           timestamps is computed, and the nearest
%                           neighbour at each timepoint in the new vector
%                           is located in the original data and copied to a
%                           resampled buffer. 
%
%   timebuffer          -   Used for any operations that need timestamp
%                           data, such as resample. 
%
%   removegaps          -   When stimulus presentation is paused, or
%                           during calibration, or is split sessions have
%                           been combined, there can be gaps in the gaze
%                           data. This is not a problem in itself, but
%                           causes issues for resampling, where the nearest
%                           temporal neighbour for a number of samples at
%                           the new sampling rate is one sample in the
%                           original buffer. This in turn leads to runs of
%                           identical timestamps (i.e. inter-sample time
%                           deltas of 0). This can cause problems when data
%                           are smoothed, filtered or interpolated, since
%                           an assumption of such methods is that
%                           timestamps are positively monotonically
%                           increasing. Specifying the removegaps option,
%                           in concert with the resampleFs operation, will
%                           remove these gaps, leading to correct,
%                           monotonic, inter-sample deltas). removegaps
%                           defaults to TRUE. 
%
%   interpolate             Linearly interpolate gaps of <200ms. See
%                           etInterpBuffer for more details/options. This
%                           option requires a timebuffer argument. 
%
%   interpolateCritMs       Criterion for interpolation, in ms. Gaps of
%                           missing data < criterion will be interpolated.
%                           Default is 200ms. 
%
%   removesaccades          A basic inter-sample velocity filter. Samples
%                           with a velocity > 50deg/s are classified as
%                           saccades. Gaze, pupil and head position data
%                           for these samples are replaced with NaNs. Note
%                           that in order to calculate metrics in degress
%                           of visual arc, you must specify 'screenwidth'
%                           and 'screenheigh' parameters (see below). Since
%                           this processing stage also needs to calculate
%                           degrees per second, you must also pass a
%                           timeBuffer parameter. 
%                           This parameter assumes that the gaze data in mb are
%                           already in degrees. If this is not the case,
%                           also call the 'converttodegrees' parameter. 
%
%   converttodegrees        Convert normalized gaze coords (in the range
%                           0-1) to degrees of visual angle. To use this
%                           option, you must specify 'screenwidth' and
%                           'screenheight' parameters (see below).
%
%   tempconverttodegrees    If you wish to perform operations that require
%                           data in degrees of visual angle (such as
%                           removing saccades according to a velocity
%                           threshold specified in degrees), but you wish
%                           to continue analysing the data in a different
%                           coordinate system (e.g. normalised), set this
%                           option to true. This will convert to degrees
%                           for processing, but return results in the
%                           original coordinate system. 
%
%   screenWidth             Screen width in cm.
%
%   screenHeight            Screen height in cm. 
%
%[mb, info] = ETPREPROCESS returns the optional output argument info. This
%contains the following fields:
% 
%   propValL            -   Proportion of valid (non-missing) samples for
%                           left eye
%
%   propValR            -   Proportion of valid samples for the right eye
%
%   propVal             -   Proportion of valid samples for either the left
%                           or the right eye
%
%
% Version 1.1 20180123


parser          = inputParser;
addParameter(   parser, 'averageeyes',          false,      @islogical      )
addParameter(   parser, 'removemissing',        false,      @islogical      )
addParameter(   parser, 'removeoffscreen',      false,      @islogical      )
addParameter(   parser, 'binocularonly',        false,      @islogical      )
addParameter(   parser, 'resampleFs',           [],         @isnumeric      )
addParameter(   parser, 'timebuffer',           [],         @isnumeric      )
addParameter(   parser, 'removegaps',           false,      @islogical      )
addParameter(   parser, 'interpolate',          false,      @islogical      )
addParameter(   parser, 'interpolateCritMs',    200,        @isnumeric      )
addParameter(   parser, 'removesaccades',       false,      @islogical      )
addParameter(   parser, 'converttodegrees',     false,      @islogical      )
addParameter(   parser, 'tempconverttodegrees', false,      @islogical      )
addParameter(   parser, 'screenwidth',          false,      @isnumeric      )
addParameter(   parser, 'screenheight',         false,      @isnumeric      )

parse(                  parser, varargin{:})
averageeyes             = parser.Results.averageeyes;
removemissing           = parser.Results.removemissing;
removeoffscreen         = parser.Results.removeoffscreen;
binocularonly           = parser.Results.binocularonly;
resampleFs              = parser.Results.resampleFs;
tb                      = parser.Results.timebuffer;
removeGaps              = parser.Results.removegaps;
interpolate             = parser.Results.interpolate;
interCritMs             = parser.Results.interpolateCritMs;
removesaccades          = parser.Results.removesaccades;
converttodegrees        = parser.Results.converttodegrees;
tempconverttodegrees    = parser.Results.tempconverttodegrees;
screenwidth             = parser.Results.screenwidth;
screenheight            = parser.Results.screenheight;

% check input args 
if isempty(resampleFs) && removeGaps
    error('removegaps can only be specified in conjunction with resampleFs.')
end

if (removesaccades || converttodegrees || tempconverttodegrees) &&...
        (isempty(screenwidth) || isempty(screenheight))
    error('If ''removesaccades'', ''converttodegrees'', or ''tempconverttodegrees'' are specified, you must also specify ''screenwidth'' and ''screenheight''.')
end

if removesaccades && isempty(tb)
    error('The ''removesaccades'' option requires a ''timeBuffer'' argument.')
end

if converttodegrees && tempconverttodegrees
    error('Cannot specify ''converttodegrees'' and ''tempconverttodegrees'' at the same time.')
end

if interpolate && isempty(tb)
    error('You must provide a timebuffer in order to interpolate.')
end

% calculate
[missL, missR, osL, osR, bino] = calculate(mb);
    
% operations
if removemissing
    mb(missL, 1:12)         = nan;
    mb(missL, 13)           = 4;
    mb(missR, 14:25)        = nan;
    mb(missR, 26)           = 4;
    [missL, missR, osL, osR, bino] = calculate(mb);
end

if removeoffscreen
    mb(osL, 7:8)            = nan;
    mb(osL, 13)             = 4;
    mb(osR, 20:21)          = nan;
    mb(osR, 26)             = 4;
    [missL, missR, osL, osR, bino] = calculate(mb);
end

if binocularonly
    mb(~bino, 1:12)         = nan;
    mb(~bino, 13)           = 4;
    mb(~bino, 14:25)        = nan;
    mb(~bino, 26)           = 4;
    [missL, missR, osL, osR, bino] = calculate(mb);
end

if averageeyes
    % average bino (not strict bino)
    notStrictBino                   = ~missL & ~missR;
    mb(notStrictBino, 7)            = mean(mb(notStrictBino, [7, 20]), 2);   % left X
    mb(notStrictBino, 8)            = mean(mb(notStrictBino, [8, 21]), 2);   % left y
    mb(notStrictBino, 20)           = mb(notStrictBino, 7);                  % right x
    mb(notStrictBino, 21)           = mb(notStrictBino, 8);                  % right y
    mb(notStrictBino, 1)            = mean(mb(notStrictBino, [1, 14]), 2);   % headpos left x
    mb(notStrictBino, 2)            = mean(mb(notStrictBino, [2, 15]), 2);   % headpos left y
    mb(notStrictBino, 3)            = mean(mb(notStrictBino, [3, 16]), 2);   % headpos left z
    mb(notStrictBino, 14)           = mb(notStrictBino, 1);                  % headpos right x
    mb(notStrictBino, 15)           = mb(notStrictBino, 2);                  % headpos right y
    mb(notStrictBino, 16)           = mb(notStrictBino, 3);                  % headpos right z   
    % copy mono to missing eye
    mb(missL & ~missR, [7, 8])      = mb(missL & ~missR, [20, 21]);     % left x, y
    mb(missL & ~missR, [1, 2, 3])   = mb(missL & ~missR, [14, 15, 16]); % left headpos x, y, z
    mb(missL & ~missR, 13)          = 0;                                % left val
    mb(missR & ~missL, [20, 21])    = mb(missR & ~missL, [7, 8]);       % right x, y
    mb(missR & ~missL, [14, 15, 16])= mb(missR & ~missL, [1, 2, 3]);    % right headpos x, y, z
    mb(missR & ~missL, 26)          = 0;                                % right val
    % recalc
    [missL, missR, osL, osR, bino]  = calculate(mb);
end

if interpolate
    cfg.dontfilter = true;
    cfg.maxms = interCritMs;
    cfg.mainbuffer = mb;
    cfg.timebuffer = tb;
    cfg.makeplot = true;
    [mb, info.interpFlags] = etInterpBuffer(cfg);
end

if converttodegrees || tempconverttodegrees
    % if only temporary, make a copy of the data in the original coordinate
    % system
    if tempconverttodegrees, mb_orig = mb; end
    % calculate median distance from screen
    mb_dis = etPreprocess(mb, 'removemissing', true,...
        'averageeyes', true, 'binocularonly', true);
    head_dis = etReturnGaze(mb_dis, 'analyticssdk', 'ldis');
    medDist = nanmedian(head_dis) / 10;
    % convert to degrees
    [mb(:, 7), mb(:, 8)] = norm2deg(mb(:, 7), mb(:, 8), screenwidth,...
        screenheight, medDist);
    [mb(:, 20), mb(:, 21)] = norm2deg(mb(:, 20), mb(:, 21), screenwidth,...
        screenheight, medDist);
end

if removesaccades
    isSac = etClassifySaccade(mb, tb);
    mb(isSac, 1:12)         = nan;
    mb(isSac, 13)           = 4;
    mb(isSac, 14:25)        = nan;
    mb(isSac, 26)           = 4;
    % if only temporarily converting to degrees, for non-saccades, replace
    % the data with the copy taken in the original coordinate system
    if tempconverttodegrees
        mb(~isSac, :) = mb_orig(~isSac, :);
    end
    % recalc
    [missL, missR, osL, osR, bino] = calculate(mb);
end

if ~isempty(resampleFs)
    
    % check input args
    if isempty(tb)
        error('To resample, you must pass a timebuffer.')
    end
    
    % get current sampling rate
    sr = etDetermineSampleRate(tb);
    if sr == resampleFs
        warning('Desired sampling rate matched existing sampling rate.')
    end
    
    % upsampling is not supported right now - error if requested
    if resampleFs > sr
        error('Upsampling not implemented!')
    end
    
    % remove samples with zero inter-sample delta
    notMonotonic = diff(tb(:, 1)) <= 0;
    mb(notMonotonic, :) = [];
    tb(notMonotonic, :) = [];
    if any(notMonotonic)
        warning('Some (%d - %.1f%%) inter-sample timestamp deltas were not monotonically increasing and have been removed.',...
            sum(notMonotonic), (sum(notMonotonic) / length(notMonotonic)) * 100)
    end
    
    % convert timebuffer to secs
    t = etTimeBuffer2Secs(tb);
    
    % calculate inter-sample delta for requested sample rate
    isd_r = 1 / resampleFs;
    
    % make time vector for resampled data, in secs
    t_rs = 0:isd_r:t(end);
    idx = arrayfun(@(x) find(t >= x, 1, 'first'), t_rs);
    
    % median filter
    medFiltSpan = round(sr / resampleFs);
    mb(:, [1:12, 14:25]) = medfilt1(mb(:, [1:12, 14:25]), medFiltSpan);
    
    % make resampled buffers
    mb = mb(idx, :);
    tb = tb(idx, :);
    
    % detect and remove gaps
    if removeGaps
        notMono = diff(idx) <= 0;
        if all(notMono)
            % check that we aren't removing all samples
            error('After resampling, all timestamps are not monotonic.')
        elseif any(notMono)
            % remove non-monotonic samples (i.e. gaps)
            mb(notMono, :) = [];
            tb(notMono, :) = [];
        end
    end
    
    % recalc
    [missL, missR, osL, osR, bino]  = calculate(mb);
       
end

numSamps = size(mb, 1);
info.missL                  = missL;
info.missR                  = missR;
info.propValL               = sum(~missL) / numSamps;
info.propValR               = sum(~missR) / numSamps;
info.propVal                = sum(~missL | ~missR) / numSamps;
info.onScreenL              = ~osL;
info.onScreenR              = ~osR;
info.propOffScreenL         = sum(~osL) / numSamps;
info.propOffScreenR         = sum(~osR) / numSamps;
end

function [missL, missR, osL, osR, bino] = calculate(mb)
    % missing data
    missL                       = mb(:, 13) == 4;
    missR                       = mb(:, 26) == 4;

    % offscreen data
    osL                         = any(mb(:, 7:8) > 1 | mb(:, 7:8) < 0, 2);
    osR                         = any(mb(:, 20:21) > 1 | mb(:, 20:21) < 0, 2);

    % binocular data
    bino                        = mb(:, 13) == 0 & mb(:, 26) == 0;
end

% TODO
%   - Support upsampling
%   - Consider more sophisticated resampling methods
