classdef etGazeDataBino < etGazeData
        
    methods
    
        function obj = Import(obj, lx, ly, rx, ry, time, missingLeft,...
                missingRight, absent, timestamps)
        % import binocular gaze data  
        
        % check input args
        
            if ~exist('lx', 'var') || ~exist('ly', 'var') ||...
                    ~exist('rx', 'var') || ~exist('ry', 'var') ||...
                    ~exist('time', 'var')
                error('Must provide, at minimum, lx, ly, rx, ry and time.')
            end
            if ~exist('missingLeft', 'var'), missingLeft = []; end
            if ~exist('missingRight', 'var'), missingRight = []; end
            
            % ensure gaze is passed as either a vector (1 subject) or a
            % matrix (multiple subs)
            val_lx      = isvector(lx) || ismatrix(lx);
            val_ly      = isvector(ly) || ismatrix(ly);
            val_rx      = isvector(rx) || ismatrix(rx);
            val_ry      = isvector(ry) || ismatrix(ry);
            val_time    = isvector(time);
            val_missL   = isempty(missingLeft) || isvector(missingLeft) || ...
                ismatrix(missingLeft);
            val_missR   = isempty(missingRight) || isvector(missingRight) || ...
                ismatrix(missingRight);
            val_format  = all([val_lx, val_ly, val_rx, val_ry, val_time,...
                            val_missL, val_missR]);
                        
            if ~val_format
                error('lx, ly, rx, ry and time must be vectors of samples, or matrices in the form of [samples, subs].')
            end                        
            
            % ensure gaze is numeric 
            val_lx      = isnumeric(lx);
            val_ly      = isnumeric(ly);
            val_rx      = isnumeric(rx);
            val_ry      = isnumeric(ry);
            val_time    = isnumeric(time);
            val_missL   = isnumeric(missingLeft) || islogical(missingLeft);
            val_missR   = isnumeric(missingRight) || islogical(missingRight);
            val_number  = all([val_lx, val_ly, val_rx, val_ry, val_time,...
                            val_missL, val_missR]);
                        
            if ~val_number
                error('lx, ly, rx, ry and time must be numeric.')
            end    
            
            % if missing is not specified, try to build from NaNs in the x,
            % y data
            if isempty(missingLeft)
                missingLeft = false(size(lx));
                if any(localInvalidGazeCoordinates(lx, ly), 'all')
                    fprintf('Missing data from LEFT EYE derived from NaNs in gaze coords.\n')
                end
            end
            if isempty(missingRight)
                missingRight = false(size(rx));
                if any(localInvalidGazeCoordinates(rx, ry), 'all')
                    fprintf('Missing data from RIGHT EYE derived from NaNs in gaze coords.\n')
                end
            end
            missingLeft = logical(missingLeft) | ...
                localInvalidGazeCoordinates(lx, ly);
            missingRight = logical(missingRight) | ...
                localInvalidGazeCoordinates(rx, ry);
            
            % if absent is not specified, assume all samples are present
            if ~exist('absent', 'var') || isempty(absent)
                absent = false(size(lx));
%                 fprintf('Absent argument not supplied, setting .Absent to false for all samples.\n')
            end
            
            % ensure that sizes match
            val_size = isequal(size(lx), size(ly), size(rx), size(ry),...
                size(missingLeft), size(missingRight), size(absent));
            val_size = val_size & isequal(size(lx, 1), size(time, 1));
                        
            if ~val_size
                error('Number of samples in gaze data (lx, ly, rx, ry, missing, absent) and time must match.')
            end
            
            % if timestamps are not specified, use time
            if ~exist('timestamps', 'var') || isempty(timestamps)
                timestamps = time;
            end

        % directly store left/right x, y, missing, absent, time
        
            lxForAverage = double(lx);
            lyForAverage = double(ly);
            rxForAverage = double(rx);
            ryForAverage = double(ry);
            lxForAverage(missingLeft) = nan;
            lyForAverage(missingLeft) = nan;
            rxForAverage(missingRight) = nan;
            ryForAverage(missingRight) = nan;
            if obj.PreserveInvalidCoordinates
                obj.LeftX = lx;
                obj.LeftY = ly;
                obj.RightX = rx;
                obj.RightY = ry;
            else
                obj.LeftX = lxForAverage;
                obj.LeftY = lyForAverage;
                obj.RightX = rxForAverage;
                obj.RightY = ryForAverage;
            end
            obj.Time = time;
            obj.Timestamp = timestamps;
            obj.LeftMissing = missingLeft;
            obj.RightMissing = missingRight;
            obj.Absent = absent;

        % average eyes to form x, y
            
            obj.X = nanmean(cat(3, lxForAverage, rxForAverage), 3);
            obj.Y = nanmean(cat(3, lyForAverage, ryForAverage), 3);
            
            % set .Missing to be samples where both left and right are
            % missing
            obj.prMissing = obj.LeftMissing & obj.RightMissing;
            
        end
        
        function obj = ImportPupil(obj, pl, pr, leftMissing, rightMissing)
        % Import left/right pupil data and their validity masks.
        %
        % Existing callers may supply one shared missing-data mask:
        %   ImportPupil(pl, pr, pupilMissing)
        %
        % New callers may preserve independent per-eye validity:
        %   ImportPupil(pl, pr, leftPupilMissing, rightPupilMissing)

            % only import pupil if gaze data already imported
            if isempty(obj)
                error('Must import gaze data before importing pupil data.')
            end

            % If no pupil masks are supplied, use the corresponding gaze
            % masks. A single supplied mask retains the historic API and is
            % silently applied to both eyes.
            if ~exist('leftMissing', 'var') || isempty(leftMissing)
                leftMissing = obj.LeftMissing;
                rightMissing = obj.RightMissing;
            elseif ~exist('rightMissing', 'var') || isempty(rightMissing)
                rightMissing = leftMissing;
            end

            expectedSize = [obj.NumSamples, obj.NumSubjects];
            if ~isnumeric(pl) || ~isnumeric(pr) || ...
                    ~isequal(size(pl), expectedSize) || ~isequal(size(pr), expectedSize)
                error(['Left and right pupil data must be numeric arrays of ', ...
                    'size [%d x %d].'], expectedSize(1), expectedSize(2))
            end

            leftMissing = localPupilMissingMask(leftMissing, expectedSize, 'left');
            rightMissing = localPupilMissingMask(rightMissing, expectedSize, 'right');

            % A NaN pupil diameter is missing regardless of an external
            % validity flag.
            leftMissing = leftMissing | isnan(pl);
            rightMissing = rightMissing | isnan(pr);

            obj.LeftPupil = pl;
            obj.RightPupil = pr;
            obj.LeftPupilMissing = leftMissing;
            obj.RightPupilMissing = rightMissing;
            obj.PupilMissing = leftMissing & rightMissing;

            % Average only the valid eyes. Keep the per-eye values intact so
            % that validity and raw diameters can round-trip independently.
            plForAverage = pl;
            prForAverage = pr;
            plForAverage(leftMissing) = nan;
            prForAverage(rightMissing) = nan;
            obj.Pupil = nanmean(cat(3, plForAverage, prForAverage), 3);

        end

        function buffer = ExportTaskEngine2(obj)

            buffer = nan(obj.NumSamples, 33);

            buffer(:, 1) = obj.Timestamp;
            buffer(:, 2) = obj.LeftX;
            buffer(:, 3) = obj.LeftY;
            buffer(:, 4) = ~obj.LeftMissing;
            if ~isempty(obj.LeftPupil)
                buffer(:, 5) = obj.LeftPupil;
                buffer(:, 6) = ~localStoredPupilMissing( ...
                    obj.LeftPupilMissing, localStoredPupilMissing( ...
                    obj.PupilMissing, obj.LeftMissing));
            else
                buffer(:, 6) = false(obj.NumSamples, 1);
            end
            buffer(:, 16) = false(obj.NumSamples, 1);
 
            buffer(:, 17) = obj.RightX;
            buffer(:, 18) = obj.RightY;
            buffer(:, 19) = ~obj.RightMissing;
            if ~isempty(obj.RightPupil)
                buffer(:, 20) = obj.RightPupil;
                buffer(:, 21) = ~localStoredPupilMissing( ...
                    obj.RightPupilMissing, localStoredPupilMissing( ...
                    obj.PupilMissing, obj.RightMissing));
            else
                buffer(:, 21) = false(obj.NumSamples, 1);
            end            
            buffer(:, 31) = false(obj.NumSamples, 1);
            buffer(:, 32) = obj.Timestamp;
            buffer(:, 33) = obj.Timestamp;

        end
        
        function [mb, tb, eb] = ExportTobiiAnalytics(obj)
            buffer_te2 = obj.ExportTaskEngine2;
            [mb, tb] = teConvertGaze(buffer_te2, [], 'taskengine2', 'tobiianalytics');
            if istable(obj.Events) && ~isempty(obj.Events) &&...
                    ismember('timestamp', obj.Events.Properties.VariableNames) &&...
                    ismember('data', obj.Events.Properties.VariableNames)
                eb = [num2cell([obj.Events.timestamp, obj.Events.timestamp]), obj.Events.data];
            else
                eb = [];
            end
        end
        
        function val = horzcat(~, varargin)
            
            numSubs = length(varargin);
            lens = cellfun(@(x) x.NumSamples, varargin);
            maxLen = max(lens);
            
            lx = nan(maxLen, numSubs);
            ly = nan(maxLen, numSubs);
            rx = nan(maxLen, numSubs);
            ry = nan(maxLen, numSubs);
            lm = false(maxLen, numSubs);
            rm = false(maxLen, numSubs);
            lp = nan(maxLen, numSubs);
            rp = nan(maxLen, numSubs);
            lpm = true(maxLen, numSubs);
            rpm = true(maxLen, numSubs);
            absent = true(maxLen, numSubs);
            for s = 1:numSubs
                
                s2 = lens(s);
                absent(1:s2, s) = false;
                lx(1:s2, s) = varargin{s}.LeftX;
                ly(1:s2, s) = varargin{s}.LeftY;
                rx(1:s2, s) = varargin{s}.RightX;
                ry(1:s2, s) = varargin{s}.RightY;
                lm(1:s2, s) = varargin{s}.LeftMissing;
                rm(1:s2, s) = varargin{s}.RightMissing;
                lp(1:s2, s) = varargin{s}.LeftPupil;
                rp(1:s2, s) = varargin{s}.RightPupil;
                lpm(1:s2, s) = localStoredPupilMissing( ...
                    varargin{s}.LeftPupilMissing, varargin{s}.PupilMissing);
                rpm(1:s2, s) = localStoredPupilMissing( ...
                    varargin{s}.RightPupilMissing, varargin{s}.PupilMissing);
                
            end
            
            idx_longest = find(lens == maxLen, 1);
            t = varargin{idx_longest}.Time;
            
            val = etGazeDataBino;
            val.Import(lx, ly, rx, ry, t, lm, rm, absent);
            if ~all(isnan(lp(:)))
                val.ImportPupil(lp, rp, lpm, rpm)
            end

        end
             
    end
        
end

function mask = localPupilMissingMask(mask, expectedSize, eyeName)

if isscalar(mask)
    mask = repmat(logical(mask), expectedSize);
elseif ~isequal(size(mask), expectedSize)
    error('%s pupil missing mask must be scalar or size [%d x %d].', ...
        eyeName, expectedSize(1), expectedSize(2))
elseif ~(islogical(mask) || isnumeric(mask))
    error('%s pupil missing mask must be logical or numeric.', eyeName)
else
    mask = logical(mask);
end

end


function mask = localStoredPupilMissing(storedMask, fallbackMask)

if isempty(storedMask)
    mask = logical(fallbackMask);
else
    mask = logical(storedMask);
end

end


function missing = localInvalidGazeCoordinates(x, y)

missing = ~isfinite(x) | ~isfinite(y) | ...
    x < 0 | x > 1 | y < 0 | y > 1;

end
