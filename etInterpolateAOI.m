function [in, postInterpMissing, postInterpX, postInterpY] =...
    etInterpolateAOI(in, gaze, maxS)
    
% validate inputs

    % if maximum gap length to fill-in is not passed, assume not upper
    % limit and set value to inf
    if ~exist('maxS', 'var') || isempty(maxS)
        maxS = inf;
    end
    
    % absent means data was not collected (as opposed to being missing).
    % Absent data cannot therefore be missing. Warn on this condition
    if any(gaze.Absent & gaze.Missing)
        gaze.Missing(gaze.Absent) = false;
        warning('Some samples are marked as missing AND absent. Absent will take priority over missing.')
    end
    
%     % data in the AOI (in the variable 'in') should not be absent or
%     % missing. If so, 'in' takes priority
%     if any(in & gaze.Absent)
%         gaze.Absent(in) = false;
%         warning('Some samples are marked as in the AOI AND absent. In the AOI will take priority over absent.')
%     end
%     if any(in & gaze.Missing)
%         gaze.Missing(in) = false;
%         warning('Some samples are marked as in the AOI AND missing. In the AOI will take priority over missing.')
%     end
    
% setup 

    % make a copy of PropValid from the gaze data. We'll update this with
    % post-interp missing vector (since we'll be filling in some missing
    % values in the AOI vector, it makes sense to at least record this)
    postInterpMissing = gaze.Missing;
    postInterpX = gaze.X;
    postInterpY = gaze.Y;
    
    % if all in, or all out, return
    if all(in(:)) || ~any(in(:))
        return
    end

    numAOIs = size(in, 3);   
    numSubs = size(in, 2);
    
    % notInAny means that gaze was not in any AOI and was not missing. Only
    % those samples can be interpolated
    notInAny = ~gaze.Absent & (all(~in, 3) | gaze.Missing);
    
% process each subject

    for s = 1:numSubs

        % find runs of missing or out-of-AOI samples
        ctm = findcontig2(gaze.Missing(:, s) | notInAny(:, s), true);
        if isempty(ctm), continue, end

        % convert to secs
        [ctm_time, ctm]         = contig2time(ctm, gaze.Time);

        % remove gaps longer than criterion 
        tooLong                 = ctm_time(:, 3) > maxS;
        ctm(tooLong, :)         = [];
        ctm_time(tooLong, :)    = [];
        dur                     = ctm_time(:, 3);

        % find samples on either side of the edges of missing data
        e1                      = ctm(:, 1) - 1;
        e2                      = ctm(:, 2) + 1;

        % remove out of bounds 
        outOfBounds             = e1 == 0 | e2 > size(in, 1);
        e1(outOfBounds)         = [];
        e2(outOfBounds)         = [];
        ctm(outOfBounds, :)     = [];
        dur(outOfBounds)        = [];

        % check each edge and flag whether a) gaze was in an AOI at both edges,
        % and b) gaze was in the SAME AOI at both edges
        val = false(length(dur), numAOIs);
        for e = 1:length(e1)

            % get state of all AOIs at edge samples
            check1 = in(e1(e), s, :);
            check2 = in(e2(e), s, :);

            % check state is valid 
            val(e, :) = sum([check1; check2], 1) == 2;

            % fill in gaps
            for a = 1:numAOIs
                if val(e, a)
                    
                    % fill in true values to 'in' AOI score vector
                    in(e1(e):e2(e), s, a) = true;
                    
                    % ensure that interpolated samples are not marked as
                    % missing
                    postInterpMissing(e1(e):e2(e), s) = false;
                    
                    % interpolate gaze coords
                    gx = gaze.X(e1(e):e2(e), s);
                    gy = gaze.Y(e1(e):e2(e), s);
                    gt = gaze.Time(e1(e):e2(e));
                    xi = interp1([gt(1), gt(end)], [gx(1), gx(end)], gt);
                    yi = interp1([gt(1), gt(end)], [gy(1), gy(end)], gt);
                    
                    % store
                    postInterpX(e1(e):e2(e), s) = xi;
                    postInterpY(e1(e):e2(e), s) = yi;
                    
                end        
            end

        end
        
    end
    
    in = logical(in);

end