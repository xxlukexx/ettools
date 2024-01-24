function in = etAOITriggerTolerance(in, gaze, tolS)
% in = ETAOITRIGGERTOLERANCE(in, gaze, tolS) takes an existing AOI score
% vector ('in') and applies a tolerance threshold. This threshold is
% specified in seconds and will allow very short periods of scoring within
% the AOI to be removed. For example, a single noisy/badly calibrated
% sample may trigger the AOI erroneously. Applying a 100ms trigger
% tolerance threshold will remove this sample from the AOI score vector.
%
%   in      An AOI score vector
%
%   gaze    A teGaze instance containing gaze data
%
%   tolS    Trigger tolerance threshold in seconds
%
% Returns 'in', a modified AOI score vector. 
%
% This function will handle multiple AOI scores. The ''in'' variable should
% be organised as a [numSamples, numAOIs] matrix.
%
% It is recommended to use etInterpolateAOI on ''in'' first, to ensure that
% looks interrupted by missing data are filled in. This function will then
% remove only those noisy outliers that are unlikely to represent valid
% gaze to the AOI.

% check input args

    if ~exist('in', 'var') || (~islogical(in) && ~isnumeric(in))
        error('Must supply ''in'' as a logical vector or matrix.')
    end

    if ~exist('gaze', 'var') || ~isa(gaze, 'etGazeData')
        error('Must supply ''gaze'' as an etGazeData instance.')
    end
    
    if all(in(:)) || ~any(in(:))
        return
    end
    
%     if length(in) ~= length(gaze)
%         error('Size mismatch between ''in'' (%d samples) and ''gaze'' (%d samples).',...
%             length(in), length(gaze))
%     end
    
    if ~exist('tolS', 'var') || isempty(tolS)
        tolS = 0.050;
        fprintf('Default trigger tolerance threshold of %.3fs used.\n', tolS);
    end
    
% process

    % loop through each AOI
    numAOIs = size(in, 3);
    for a = 1:numAOIs
        
        numSubs = size(in, 2);
        for s = 1:numSubs
        
            if all(in(:, s, a)) || ~any(in(:, s, a))
                continue
            end

             % find contiguous runs of samples
            ct = findcontig2(in(:, s, a));
            ctt = contig2time(ct, gaze.Time);

            % find runs below trigger threshold
            idx_trig = ctt(:, 3) < tolS;

            % remove runs ABOVE threshold (we only work with those
            % below)
            ct = ct(idx_trig, :);
            numTrig = sum(idx_trig);

            % convert below-threshold runs back to logical index
            for trig = 1:numTrig

                % find sample edges of look that is being deleted on
                % account of being below trigger threshold
                s1 = ct(trig, 1);
                s2 = ct(trig, 2);

                % delete look
                in(s1:s2, s, a) = false;

            end
            
        end

    end

end