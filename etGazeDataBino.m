classdef etGazeDataBino < etGazeData
        
    methods
        
%         function obj = etGazeDataBino(varargin)
%             obj = obj@etGazeData;
%         end
    
        function obj = Import(obj, lx, ly, rx, ry, time, missingLeft,...
                missingRight, absent, timestamps)
        % import binocular gaze data  
        
        % check input args
        
            if ~exist('lx', 'var') || ~exist('ly', 'var') ||...
                    ~exist('rx', 'var') || ~exist('ry', 'var') ||...
                    ~exist('time', 'var')
                error('Must provide, at minimum, lx, ly, rx, ry and time.')
            end
            
            % ensure gaze is passed as either a vector (1 subject) or a
            % matrix (multiple subs)
            val_lx      = isvector(lx) || ismatrix(lx);
            val_ly      = isvector(ly) || ismatrix(ly);
            val_rx      = isvector(rx) || ismatrix(rx);
            val_ry      = isvector(ry) || ismatrix(ry);
            val_time    = isvector(time);
            val_missL   = isvector(missingLeft) || ismatrix(missingLeft);
            val_missR   = isvector(missingRight) || ismatrix(missingRight);
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
            if ~exist('missingLeft', 'var') || isempty(missingLeft)
                missingLeft = isnan(lx) | isnan(ly);
                if any(missingLeft)
                    fprintf('Missing data from LEFT EYE derived from NaNs in gaze coords.\n')
                end
            end
            if ~exist('missingRight', 'var') || isempty(missingRight)
                missingRight = isnan(rx) | isnan(ry);
                if any(missingRight)
                    fprintf('Missing data from RIGHT EYE derived from NaNs in gaze coords.\n')
                end
            end
            
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
        
            obj.LeftX = lx;
            obj.LeftY = ly;
            obj.RightX = rx;
            obj.RightY = ry;
            obj.Time = time;
            obj.Timestamp = timestamps;
            obj.LeftMissing = missingLeft;
            obj.RightMissing = missingRight;
            obj.Absent = absent;

        % average eyes to form x, y
            
            obj.X = nanmean(cat(3, lx, rx), 3);
            obj.Y = nanmean(cat(3, ly, ry), 3);
            
            % set .Missing to be samples where both left and right are
            % missing
            obj.prMissing = obj.LeftMissing & obj.RightMissing;
            
        end
        
        function obj = ImportPupil(obj, pl, pr, pinvalid)
            
%             warning('Beta - this does not have all the error checking it needs.')
            
            % only import pupil if gaze data already imported
            if isempty(obj)
                error('Must import gaze data before importing pupil data.')
            end
            
            % if no valid passed, using the gaze missing variables
            if ~exist('pinvalid', 'var') || isempty(pinvalid)
                pinvalid = obj.Missing;
            end
            
%             % check that all import args are same size
%             if (~isscalar(pl) || ~isscalar(pr) || ~isscalar(pinvalid)) &&...
%                     (~isvector(pl) || ~isvector(pr) || ~isvector(pinvalid))
%                 error('All pupil data must be either scalar or vectors, and all must be the same size.')
%             end
            
            % check that pupil and gaze data sizes match
            if size(pl, 1) ~= obj.NumSamples
                error('Gaze data (%d samples) / pupil data (%d samples) size mismatch.',...
                    obj.NumSamples, size(pl, 1))
            end
            
            % check that pupil and gaze data sizes match
            if size(pl, 2) ~= obj.NumSubjects
                error('Gaze data (%d subjects) / pupil data (%d subjects) size mismatch.',...
                    obj.NumSubjects, size(pl, 2))
            end
            
            obj.LeftPupil = pl;
            obj.RightPupil = pr;
            obj.PupilMissing = pinvalid;
            
            % average left/right pupil
            obj.Pupil = nanmean(cat(3, pl, pr), 3);
            
        end
        
        function buffer = ExportTaskEngine2(obj)
            
            buffer = nan(obj.NumSamples, 33);
            
            buffer(:, 1) = obj.Timestamp;
            buffer(:, 2) = obj.LeftX;
            buffer(:, 3) = obj.LeftY;
            buffer(:, 4) = ~obj.LeftMissing;
            if ~isempty(obj.LeftPupil)
                buffer(:, 5) = obj.LeftPupil;
                buffer(:, 6) = ~obj.LeftMissing;
            else
                buffer(:, 6) = false(obj.NumSamples, 1);
            end
            buffer(:, 16) = false(obj.NumSamples, 1);
 
            buffer(:, 17) = obj.RightX;
            buffer(:, 18) = obj.RightY;
            buffer(:, 19) = ~obj.RightMissing;
            if ~isempty(obj.RightPupil)
                buffer(:, 20) = obj.RightPupil;
                buffer(:, 21) = ~obj.RightMissing;
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
            pm = false(maxLen, numSubs);
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
                pm(1:s2, s) = varargin{s}.PupilMissing;
                
            end
            
            idx_longest = find(lens == maxLen, 1);
            t = varargin{idx_longest}.Time;
            
            val = etGazeDataBino;
            val.Import(lx, ly, rx, ry, t, lm, rm, absent);
            if ~all(isnan(lp(:)))
                val.ImportPupil(lp, rp, pm)
            end

        end
             
    end
        
end
