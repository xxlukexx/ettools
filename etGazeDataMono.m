classdef etGazeDataMono < etGazeData

    methods 

        function obj = Import(obj, x, y, time, missing, absent)
        % import monocular gaze data  
        
            % check input args
            if ~exist('x', 'var') || ~exist('y', 'var') ||...
                    ~exist('time', 'var')
                error('Must provide, at minimum, x, y, time.')
            end
            if ~exist('missing', 'var'), missing = []; end
            
            if ~isvector(x) || ~isvector(y) || ~isvector(time)
                error('x, y, and time must be vectors.')
            end
            
            if ~isnumeric(x) || ~isnumeric(y) || ~isnumeric(time)
                error('x, y, and time must be numeric.')
            end
            
            if ~isequal(size(x), size(y), size(time))
                error('x, y, and time must be of the same size.')
            end
            
            if isempty(missing)
                missing = false(size(x));
            elseif ~(isnumeric(missing) || islogical(missing)) || ...
                    ~isequal(size(missing), size(x))
                error('Missing mask must be numeric or logical and match x/y size.')
            end
            missing = logical(missing) | localInvalidGazeCoordinates(x, y);

            xForAnalysis = double(x);
            yForAnalysis = double(y);
            xForAnalysis(missing) = nan;
            yForAnalysis(missing) = nan;
            if obj.PreserveInvalidCoordinates
                obj.X = x;
                obj.Y = y;
            else
                obj.X = xForAnalysis;
                obj.Y = yForAnalysis;
            end

            % left and right are copies for monocular gaze data
            obj.LeftX = obj.X;
            obj.LeftY = obj.Y;
            obj.RightX = obj.X;
            obj.RightY = obj.Y;
            obj.Time = time;
            obj.Missing = missing;
            
            % duplicate Missing for left and right eye
            obj.LeftMissing = obj.Missing;
            obj.RightMissing = obj.Missing;
            
            % if absent is not specified, assume all samples are present
            if ~exist('absent', 'var') || isempty(absent)
                obj.Absent = false(size(x));
                fprintf('Absent argument not supplied, setting .Absent to false for all samples.\n')
            else
                obj.Absent = absent;
            end
            
        end
        
        function obj = ImportPupil(obj, p, pinvalid)
            
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
            if size(p, 1) ~= obj.NumSamples
                error('Gaze data (%d samples) / pupil data (%d samples) size mismatch.',...
                    obj.NumSamples, size(p, 1))
            end

            obj.Pupil = p;
            
        end

   end
   
end


function missing = localInvalidGazeCoordinates(x, y)

missing = ~isfinite(x) | ~isfinite(y) | ...
    x < 0 | x > 1 | y < 0 | y > 1;

end
