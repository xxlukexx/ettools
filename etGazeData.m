classdef etGazeData < dynamicprops
    
    properties
        X
        Y
        LeftX
        LeftY
        RightX
        RightY
        Time
        Timestamp
        LeftMissing
        RightMissing
        Absent
        Pupil
        LeftPupil 
        RightPupil
        PupilMissing
        Events
        ScreenDimensions = [nan, nan]
        DistanceFromScreen = nan
    end
    
    properties (Dependent)
        Missing
    end
    
    properties (Dependent, SetAccess = private)
        Duration
        NumSubjects
        NumSamples
        NumValid
        LeftNumValid
        RightNumValid
        PropValid
        LeftPropValid
        RightPropValid
        TimeValid
        LeftTimeValid
        RightTimeValid
        XScreen
        YScreen
        XDeg
        YDeg
        MeanSampleRate
        HasGaze
        HasPupil
        HasEvents
        Calc
    end
    
    properties (Access = protected)        
        prMissing
        prAddedVariables
    end
    
    methods
        
        function obj = etGazeData(varargin)
        % this function is common to mono and bino gaze data (i.e. it will
        % be called when EITHER an etGazeDataMono OR etGazeDataBino is
        % instantiated). If bino data is passed to an etGazeDataMono
        % instance then the eyes will be averaged first. This provides a
        % mechanism to easily load bino data but to use the averaged eye
        % coords. Useful for noisy gaze data when left/right eyes may
        % straddle two AOIs, which causes probems when scoring. 
            
            obj.prAddedVariables = teCollection;
            
            parser = inputParser;
            parser.addParameter('mainBuffer', []);
            parser.addParameter('timeBuffer', []);
            parser.addParameter('eventBuffer', []);
            parser.addParameter('te1_preproc', []);
            parser.addParameter('te2', []);
            parser.parse(varargin{:});

            % mainbuffer and timebuffer
            mb = parser.Results.mainBuffer;
            tb = parser.Results.timeBuffer;
            eb = parser.Results.eventBuffer;
            if isempty(tb) && ~isempty(mb)
                error('If supplying a mainBuffer, must also supply a timeBuffer.')
            elseif isempty(mb) && ~isempty(tb)
                error('If supplying a timeBuffer, must also supply a mainBuffer.')
            elseif isempty(mb) && isempty(tb) && ~isempty(eb)
                error('If supplying an eventBuffer, must also supply a mainBuffer and timeBuffer.')
            elseif ~isempty(mb) && ~isempty(tb)
                % remove missing
                mb = etPreprocess(mb, 'removemissing', true);
                % preproc events
                if ~isempty(eb)
                    t_ev = double(cell2mat(eb(:, 2)) - eb{1, 2}) / 1e6;
                    tmp = etListEvents(eb);
                    lab_ev = tmp(2:end, 4);
                    obj.Events = array2table(t_ev, 'VariableNames', {'timestamp'});
                    obj.Events.data = lab_ev;
                end
                % import gaze
                lx = mb(:, 7);
                ly = mb(:, 8);
                rx = mb(:, 20);
                ry = mb(:, 21);
                missingL = isnan(lx) | lx == -1 | isnan(ly) | ly == -1;
                missingR = isnan(rx) | rx == -1 | isnan(ry) | ry == -1;
                time = etTimeBuffer2Secs(tb);
                if isa(obj, 'etGazeDataBino')
                    % import bino gaze data 
                    obj.Import(lx, ly, rx, ry, time, missingL, missingR);
                    % import pupil
                    lp = mb(:, 12);
                    rp = mb(:, 25);
                    obj.ImportPupil(lp, rp);    
                elseif isa(obj, 'etGazeDataMono')
                    % average eyes, then import
                    x = nanmean(cat(3, lx, rx), 3);
                    y = nanmean(cat(3, ly, ry), 3);
                    missing = missingL & missingR;
                    obj.Import(x, y, time, missing);
                    % import pupil
                    lp = mb(:, 12);
                    rp = mb(:, 25);
                    p = nanmean(cat(3, lp, rp), 3);
                    obj.ImportPupil(p);                    
                else
                    error('Unexpected subclass format -- should be etGazeDataMono or etGazeDataBino.')
                end        
                % store tobii format timestamps, converted to secs
                obj.Timestamp = tb(:, 1) / 1e6;

            end
            
            % te1 preproc import
            if ~isempty(parser.Results.te1_preproc)
                data = parser.Results.te1_preproc;
                % check data format
                val_format =...
                    isstruct(data) &&...
                    hasField(data, 'ParticipantID');
                if val_format 
                    if hasField(data, 'Segmentation')
                        data = renameStructField(data, 'Segmentation', 'Segments');
                    end
                    val_format = val_format && isstruct(data.Segments) &&...
                        hasField(data.Segments, 'Segment') &&...
                        hasField(data.Segments(1).Segment, 'MainBuffer') &&...
                        hasField(data.Segments(1).Segment, 'TimeBuffer') &&...
                        hasField(data.Segments(1).Segment, 'EventBuffer'); 
                end
                if ~val_format
                    error('When using the ''te1_preproc'' option, the second argument must be a struct.')
                end
                % loop through segments
                tmp = etGazeDataBino;
                for s = 1:length(data.Segments)
                    mb = data.Segments(s).MainBuffer;
                    tb = data.Segments(s).TimeBuffer;
                    eb = data.Segments(s).EventBuffer;
                    tmp(s) = etGazeDataBino(...
                        'mainBuffer', mb,...
                        'timeBuffer', tb,...
                        'eventBuffer', eb);
                end
                obj = tmp;
            end
                    
            
            % te2 gaze import
            if ~isempty(parser.Results.te2)
                buffer = parser.Results.te2;
                if isempty(buffer)
                    error('Empty gaze buffer.')
                else
                    lx = buffer(:, 2);
                    ly = buffer(:, 3);
                    rx = buffer(:, 17);
                    ry = buffer(:, 18);
                    missingL = ~buffer(:, 4);
                    missingR = ~buffer(:, 19);
                    time = buffer(:, 1) - buffer(1, 1);
                    timestamps = buffer(:, 1);
                    obj.Import(lx, ly, rx, ry, time, missingL, missingR, [], timestamps);
                end
            end            
            
        end
        
        function val = FilterOneSample(obj, s)
            
            if ~isempty(obj.prAddedVariables)
                warning('AddedVariables are not currently carried over when this method is called. This needs adding!')
            end
            
            if isa(obj, 'etGazeDataMono')
                
                val = etGazeDataMono;
                val.Import(obj.X(s, :), obj.Y(s, :), obj.Time(s, :),...
                    obj.Missing(s, :), obj.Absent(s, :));
                
            elseif isa(obj, 'etGazeDataBino')
                
                val = etGazeDataBino;
                val.Import(obj.LeftX(s, :), obj.LeftY(s, :), obj.RightX(s, :),...
                    obj.RightY(s, :), obj.Time(s, :), obj.LeftMissing(s, :),...
                    obj.RightMissing(s, :), obj.Absent(s, :));
                
                if obj.HasPupil
                    val.ImportPupil(obj.LeftPupil(s, :), obj.RightPupil(s, :),...
                        obj.PupilMissing(s, :));
                end
                
            end
  
        end
        
        function val = FilterOneSubject(obj, sub)
            
            if ~isempty(obj.prAddedVariables)
                warning('AddedVariables are not currently carried over when this method is called. This needs adding!')
            end
            
            if isa(obj, 'etGazeDataMono')
                
                val = etGazeDataMono;
                val.Import(obj.X(:, sub), obj.Y(:, sub), obj.Time(:, sub),...
                    obj.Missing(:, sub), obj.Absent(:, sub));
                
            elseif isa(obj, 'etGazeDataBino')
                
                val = etGazeDataBino;
                val.Import(obj.LeftX(:, sub), obj.LeftY(:, sub), obj.RightX(:, sub),...
                    obj.RightY(:, sub), obj.Time, obj.LeftMissing(:, sub),...
                    obj.RightMissing(:, sub), obj.Absent(:, sub));
                if obj.HasPupil
                    val.ImportPupil(obj.LeftPupil(:, sub), obj.RightPupil(:, sub),...
                        obj.PupilMissing(:, sub));
                end
                
            end
  
        end
        
        function val = SegmentBySample(obj, s1, s2)
        % take two sample indices (s1, s2) and spawn a new teGazeData 
        % instance containing just the gaze within those samples. s1 and s2
        % can be vectors of indices, in which case return a cell array of
        % teGaze instances
        
            if ~isnumeric(s1) || ~isnumeric(s2)
                error('Sample indices must be numeric.')
            end
            
            if ~isequal(size(s1), size(s2))
                error('s1 and s2 must be the same size.')
            end
            
            if ~isscalar(s1) && ~isvector(s1)
                error('s1 must be a scalar or vector of sample indices.')
            end
            
            if ~isscalar(s2) && ~isvector(s2)
                error('s2 must be a scalar or vector of sample indices.')
            end
            
            if any(s1 < 1) || any(s2 < 1) || any(s1 > obj.NumSamples) ||...
                    any(s2 > obj.NumSamples)
                error('Index out of bounds: s1 = %d, s2 = %d, size of gaze = %d',...
                    s1, s2, obj.NumSamples)
            end
            
            if s1 > s2
                error('s1 must be less than or equal to s2.')
            end
        
            numSegs = length(s1);
            val = cell(numSegs, 1);
            for s = 1:numSegs
                
                idx = s1(s):s2(s);
                    
                if isa(obj, 'etGazeDataMono')

                    val{s} = etGazeDataMono;
                    val{s}.Import(...
                        obj.X(idx, :),...
                        obj.Y(idx, :),...
                        obj.Time(idx, :),...
                        obj.Missing(idx, :),...
                        obj.Absent(idx, :));
                    if obj.HasPupil
                        error('Pupil data not yet supported for monocular data.')
                        val{s}.ImportPupil(...
                            obj.Pupil(idx, :),...
                            obj.PupilMissing(idx, :));
                    end

                elseif isa(obj, 'etGazeDataBino')

                    val{s} = etGazeDataBino;
                    val{s}.Import(...
                        obj.LeftX(idx, :),...
                        obj.LeftY(idx, :),...
                        obj.RightX(idx, :),...
                        obj.RightY(idx, :),...
                        obj.Time(idx, :) - obj.Time(idx(1), :),...  % zero time
                        obj.LeftMissing(idx, :),...
                        obj.RightMissing(idx, :),...
                        obj.Absent(idx, :),...
                        obj.Timestamp(idx, :));
                    if obj.HasPupil
                        val{s}.ImportPupil(...
                            obj.LeftPupil(idx, :),...
                            obj.RightPupil(idx, :),...
                            obj.PupilMissing(idx, :));
                    end
                    
                end        
                
                % segment events
                if obj.HasEvents
                    t1 = obj.Sample2Time(s1);
                    t2 = obj.Sample2Time(s2);
                    idx_ev = obj.Events.Time >= t1 & obj.Events.Time <= t2;
                    val{s}.Events = obj.Events(idx_ev, :);
                end
                
                % handle added variables
                for v = 1:obj.prAddedVariables.Count
                    
                    % read in the variable's value
                    addVar = obj.prAddedVariables.Keys{v};
                    addVal = obj.(addVar);
                    
                    % segment
                    addVal = addVal(idx, :);
                    
                    % write to new obj
                    val{s}.AddVariable(addVar, addVal)
                    
                end
                
                % store segment onset in output gaze object, both in
                % samples and in time
                addprop(val{s}, 'SegmentOnsetSamp');
                addprop(val{s}, 'SegmentOnsetTime');
                addprop(val{s}, 'SegmentOffsetSamp');
                addprop(val{s}, 'SegmentOffsetTime');
                
                val{s}.SegmentOnsetSamp = s1(s);
                val{s}.SegmentOnsetTime = obj.Time(s1(s));
                val{s}.SegmentOffsetSamp = s2(s);
                val{s}.SegmentOffsetTime = obj.Time(s2(s));
                
            end
            
            if numel(val) == 1
                val = val{1};
            end
            
        end
        
        function val = SegmentByTime(obj, t1, t2)
        % takes two timestamps (t1, t2) representing the edges of the
        % desired segment. Converts time to samples and then calls
        % SegmentBySample to actually chop up the data
        
            s1 = obj.Time2Sample(t1);
            s2 = obj.Time2Sample(t2);
            val = obj.SegmentBySample(s1, s2);
        
        end
        
        % dynamic props
        function AddVariable(obj, var, val)
        % adds a new variable to the object. This variable must be of size
        % [numSubjects, numSamples] in the existing object. The data is
        % returned by the SegmentBySample method. This allows the user to
        % add additional timeseries data and have it easily segmented along
        % with the gaze data.
        
            % check inputs
            if ~ischar(var) || ~isvarname(var)
                error('''var'' must be a valid Matlab variable name.')
            end
            if ~isequal(size(val), size(obj.X))
                error('The added variable must be of the same [numSubject, numSamples] size as the gaze object.')
            end
            
            obj.prAddedVariables(var) = addprop(obj, var);
            obj.(var) = val;
        end
        
        % plotting
        function PlotGaze(obj, h)
            
            if ~exist('h', 'var') || isempty(h)
                fig = figure;
            else
                fig = h;
            end
            
            if isa(fig, 'figure')
                whitebg(fig, 'k')
            end
            
            % plot missing data
            hold on
            set(gca, 'ydir', 'reverse')
            ctt = contig2time(findcontig2(obj.Missing, true), obj.Time);
            for i = 1:size(ctt, 1)
                x1 = ctt(i, 1);
                x2 = ctt(i, 3);
                y1 = 0.00;
                y2 = 1.00;
                rectangle('position', [x1, y1, x2, y2], 'FaceColor', [1, 0, 0, .1], 'EdgeColor', [1, 0, 0, .9])
                rectangle('position', [x1, 0.98, x2, 0.02], 'FaceColor', [1, 0, 0, .9], 'EdgeColor', [1, 0, 0, .9])
            end
            
            % plot x, y
            sc = scatter(obj.Time, obj.X);
            sc.MarkerFaceColor = sc.MarkerEdgeColor;
            sc.MarkerFaceAlpha = 0.5;
            hold on
            sc = scatter(obj.Time, obj.Y);
            sc.MarkerFaceColor = sc.MarkerEdgeColor;
            sc.MarkerFaceAlpha = 0.5;            
            xlabel('Time (s)')
            ylabel('X/Y Position (normalised coords)')
            ylim([0, 1])
            xlim([obj.Time(1), obj.Time(end)])
            set(gca, 'XGrid', 'on')
            set(gca, 'YGrid', 'on')
            
            % plot events
            tb = [0, 0, 0, 0];
            x = 0;
            tx = [];
            if ~isempty(obj.Events)
                ty = 0;
                for i = 1:size(obj.Events, 1)
                    ox = x;
                    x = obj.Events.timestamp(i);
                    lab = obj.Events.data{i};
                    if contains(lab, 'FRAME_CALC') ||...
                       contains(lab, 'NATSCENES_FRAME')
                        continue
                    end
                    line([x, x], [0, 1], 'Color', 'w')
                    if ~isempty(tx)
                        moveDown = x < tx.Extent(1) + tx.Extent(3);
                    else
                        moveDown = false;
                    end
                    tb = textBounds(lab, gca);
                    if moveDown
                        ty = ty + tb(4);
                    end
                    tx = text(x, ty, lab, 'Interpreter', 'none');
                end
            end            
            
            legend('X', 'Y')
            
%             % attemp to plot AOIs, in any found in events
%             idx_aoi = cellfun(@(x) contains(x, 'ADDAOI'), obj.Events.Label);
%             ev_aoi = obj.Events(idx_aoi, :);
%             for a = 1:numAOIs
%                 
%                 parts = strsplit(ev_aoi.Label{a}, '_');
%                 rect = [...
%                     str2double(parts{2}),...
%                     str2double(parts{3}),...
%                     str2double(parts{4}),...
%                     str2double(parts{5}),...
%                     ];
%                 
%                 
%                 
%             end

            % details
            pos = get(gca, 'InnerPosition');
            x1 = pos(1);
            w = pos(3) / 3;
            h = (1 - pos(4)) / 2;
            x2 = pos(1) + w;
            x3 = pos(1) + (w * 2);
            annotation('textbox', 'Position', [x1, 1 - h, w, h], 'String',...
                sprintf('Duration: %.1fs', obj.Duration),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none');
            annotation('textbox', 'Position', [x2, 1 - h, w, h], 'String',...
                sprintf('Valid: %.2f%%', obj.PropValid * 100),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none',...
                'HorizontalAlignment', 'center');
            annotation('textbox', 'Position', [x3, 1 - h, w, h], 'String',...
                sprintf('Mean fs: %.2fHz', obj.MeanSampleRate),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none',...
                'HorizontalAlignment', 'right');
            
        end
        
        function PlotPupil(obj, h, varargin)
            
            if ~exist('h', 'var') || isempty(h)
                fig = figure;
            else
                fig = h;
            end
            
            % handle options to plot only left/right/average, default to
            % plotting all three
            plotLeft = ismember(lower(varargin), '-plotleft');
            plotRight = ismember(lower(varargin), '-plotright');
            plotAverage = ismember(lower(varargin), '-plotaverage');
            if isempty(varargin) || all([~plotLeft, ~plotRight, ~plotAverage])
                plotLeft = true;
                plotRight = true;
                plotAverage = true;
            end
            
            if isa(h, 'figure')
                whitebg(fig, 'k')
            end
            
            % calculate ylim (to allow missing data, which is drawn first,
            % to be properly scaled)
            yl_min = min([obj.LeftPupil; obj.RightPupil; obj.Pupil]);
            yl_max = max([obj.LeftPupil; obj.RightPupil; obj.Pupil]);

            % plot missing data
            hold on
            set(gca, 'ydir', 'reverse')
            ctt = contig2time(findcontig2(obj.Missing, true), obj.Time);
            for i = 1:size(ctt, 1)
                x1 = ctt(i, 1);
                x2 = ctt(i, 3);
                y1 = yl_min;
                y2 = yl_max;
                rectangle('position', [x1, y1, x2, y2], 'FaceColor', [1, 0, 0, .1])
                rectangle('position', [x1, y1, x2, y2], 'FaceColor', [1, 0, 0, .9])
            end
            
            % plot left and right pupil
            hold on
            if plotLeft
                sc = scatter(obj.Time, obj.LeftPupil);
                sc.MarkerFaceColor = sc.MarkerEdgeColor;
                sc.MarkerFaceAlpha = 0.5;
            end
            if plotRight
                sc = scatter(obj.Time, obj.RightPupil);
                sc.MarkerFaceColor = sc.MarkerEdgeColor;
                sc.MarkerFaceAlpha = 0.5;     
            end
            if plotAverage
                sc = scatter(obj.Time, obj.Pupil);
                sc.MarkerFaceColor = sc.MarkerEdgeColor;
                sc.MarkerFaceAlpha = 0.5;                 
            end
            
            xlabel('Time (s)')
            ylabel('Pupil size (mm)')
            xlim([obj.Time(1), obj.Time(end)])
            ylim([yl_min, yl_max])
            set(gca, 'XGrid', 'on')
            set(gca, 'YGrid', 'on')
            
            % plot events
            tb = [0, 0, 0, 0];
            x = 0;
            tx = [];
            if ~isempty(obj.Events)
                ty = 0;
                for i = 1:size(obj.Events, 1)
                    ox = x;
                    x = obj.Events.Time(i);
                    lab = obj.Events.Label{i};
                    line([x, x], [0, 1], 'Color', 'w')
                    if ~isempty(tx)
                        moveDown = x < tx.Extent(1) + tx.Extent(3);
                    else
                        moveDown = false;
                    end
                    tb = textBounds(lab, gca);
                    if moveDown
                        ty = ty + tb(4);
                    end
                    tx = text(x, ty, lab, 'Interpreter', 'none');
                end
            end            
            
            legend('Left Pupil', 'Right Pupil', 'Avg Left/Right Pupil')

            % details
            pos = get(gca, 'InnerPosition');
            x1 = pos(1);
            w = pos(3) / 3;
            h = (1 - pos(4)) / 2;
            x2 = pos(1) + w;
            x3 = pos(1) + (w * 2);
            annotation('textbox', 'Position', [x1, 1 - h, w, h], 'String',...
                sprintf('Duration: %.1fs', obj.Duration),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none');
            annotation('textbox', 'Position', [x2, 1 - h, w, h], 'String',...
                sprintf('Valid: %.2f%%', obj.PropValid * 100),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none',...
                'HorizontalAlignment', 'center');
            annotation('textbox', 'Position', [x3, 1 - h, w, h], 'String',...
                sprintf('Mean fs: %.2fHz', obj.MeanSampleRate),...
                'VerticalAlignment', 'middle', 'Color', 'w', 'FontSize', 14,...
                'FontWeight', 'bold', 'LineStyle', 'none',...
                'HorizontalAlignment', 'right');            
            
        end
        
        % utils
        function s = Time2Sample(obj, t)
            if ~isscalar(t)
                error('Timestamp must be scalar.')
                % todo - make this not so
            elseif t > obj.Duration
                error('Requested timestamp %.2fs is greater than total duration (%.2fs)',...
                    t, obj.Duration)
            elseif t < 0 
                error('Requested timestamp is negative.')
            end
            s = find(obj.Time >= t, 1);
        end
        
        function t = Sample2Time(obj, s)
            if ~isscalar(s)
                error('Sample index must be scalar.')
                % todo - make this not so
            elseif s > obj.NumSamples
                error('Requested sample index %d is greater than total number of samples (%d)',...
                    s, obj.NumSamples)
            elseif s < 0 
                error('Requested timestamp is negative.')
            end
            t = obj.Time(s);
        end
        
        function cm = MakeCovarianceMatrix(obj, gridWidth, gridHeight)

            cm = nan(gridHeight, gridWidth);

        end
        
        % get/set
        function val = get.Duration(obj)
            val = max(obj.Time);
        end
        
        function val = get.NumSubjects(obj)
            val = size(obj.LeftX, 2);
        end
        
        function val = get.NumSamples(obj)
            val = size(obj.LeftX, 1);
        end
        
        function val = get.NumValid(obj)
            if isempty(obj)
                val = [];
            else
                val = sum(~obj.Missing & ~obj.Absent, 1);
            end
        end
        
        function val = get.LeftNumValid(obj)
            if isempty(obj)
                val = [];
            else
                val = sum(~obj.LeftMissing & ~obj.Absent, 1);
            end
        end  
        
        function val = get.RightNumValid(obj)
            if isempty(obj)
                val = [];
            else
                val = sum(~obj.RightMissing & ~obj.Absent, 1);
            end
        end          
        
        function val = get.PropValid(obj)
            if isempty(obj)
                val = [];
            else
                valid = sum(~obj.Missing & ~obj.Absent, 1);
                present = sum(~obj.Absent, 1);
                val = valid ./ present;
%                 val = sum(~obj.Missing, 1) ./ sum(~obj.Absent, 1);
            end
        end
        
        function val = get.LeftPropValid(obj)
            if isempty(obj)
                val = [];
            else
                valid = sum(~obj.LeftMissing & ~obj.Absent, 1);
                present = sum(~obj.Absent, 1);
                val = valid ./ present;                
%                 val = sum(~obj.LeftMissing, 1) ./ sum(~obj.Absent, 1);
            end
        end        
        
        function val = get.RightPropValid(obj)
            if isempty(obj)
                val = [];
            else
                valid = sum(~obj.RightMissing & ~obj.Absent, 1);
                present = sum(~obj.Absent, 1);
                val = valid ./ present;                
%                 val = sum(~obj.RightMissing, 1) ./ sum(~obj.Absent, 1);
            end
        end        
        
        function val = get.TimeValid(obj)
            if isempty(obj)
                val = [];
            elseif obj.NumSamples == 1
                val = nan;
            else
                delta = [0; diff(obj.Time)];
                val = sum(delta(~obj.Missing & ~obj.Absent));
            end
        end        
        
        function val = get.LeftTimeValid(obj)
            if isempty(obj)
                val = [];
            else
                delta = [0; diff(obj.Time)];
                val = sum(delta(~obj.LeftMissing & ~obj.Absent));
            end
        end   
        
        function val = get.RightTimeValid(obj)
            if isempty(obj)
                val = [];
            else
                delta = [0; diff(obj.Time)];
                val = sum(delta(~obj.RightMissing & ~obj.Absent));
            end
        end  
        
        function val = get.Missing(obj)
            val = obj.prMissing;
        end
        
        function set.Missing(obj, val)
            obj.prMissing = val;
            % if the overall (average both eyes) .Missing property is set,
            % apply this to both .LeftMissing and .RightMissing
            obj.LeftMissing = val;
            obj.RightMissing = val;
        end
        
        function val = get.XScreen(obj)
            val = obj.X .* obj.ScreenDimensions(1);
        end
        
        function val = get.YScreen(obj)
            val = obj.Y .* obj.ScreenDimensions(2);
        end
        
        function val = get.XDeg(obj)
            [val, ~] = norm2deg(obj.X, obj.Y, obj.ScreenDimensions(1),...
                obj.ScreenDimensions(2), obj.DistanceFromScreen);
        end
        
        function val = get.YDeg(obj)
            [~, val] = norm2deg(obj.X, obj.Y, obj.ScreenDimensions(1),...
                obj.ScreenDimensions(2), obj.DistanceFromScreen);
        end
        
        function val = get.MeanSampleRate(obj)
            val = etDetermineSampleRate(obj.Time * 1e6);
        end
        
        function val = get.HasGaze(obj)
            val = ~isempty(obj);
        end
        
        function val = get.HasPupil(obj)
            val = ~isempty(obj.Pupil);
        end
        
        function val = get.HasEvents(obj)
            val = ~isempty(obj.Events);
        end
        
        function val = get.Calc(obj)
            val = struct;
            
            % vel
            ldx = [nan; diff(obj.LeftX)];
            ldy = [nan; diff(obj.LeftY)];
            rdx = [nan; diff(obj.RightX)];
            rdy = [nan; diff(obj.RightY)];
            vl = sqrt((ldx .^ 2) + (ldy .^ 2));
            vr = sqrt((rdx .^ 2) + (rdy .^ 2));   
            
            % acc
            al = [nan; diff(vl)];
            ar = [nan; diff(vr)];
            
            % direction
            dl = 180 + atan2d(ldy, ldx);
            dr = 180 + atan2d(rdy, rdx);

            val.VelocityLeft = vl;
            val.VelocityRight = vr;
            val.VelocityLeftNormalised = vl ./ max(vl);
            val.VelocityRightNormalised = vr ./ max(vr);            
            val.AccelerationLeft = al;
            val.AccelerationRight = ar;
            val.DirectionLeft = dl;
            val.DirectionRight = dr;
            
        end

        
        % overridden methods
        function m = double(obj)
            props = properties(obj);
            m = nan(length(obj), length(props));
            for p = 1:length(props)
                m(:, p) = obj.(props{p});
            end
        end
        
        function tab = table(obj)
            if isempty(obj)
                tab = [];
            else
                tab = array2table(double(obj), 'VariableNames',...
                    properties(obj));
            end
        end
        
        function h = uitable(obj, varargin)
            h = table2uitable(table(obj), varargin{:});
        end
        
        function h = plot(obj)
            obj.PlotGaze
        end
        
        function scatter(obj)
            scatter(obj.LeftX, obj.LeftY, 5, 'g');
            hold on
            scatter(obj.RightX, obj.RightY, 5, 'b');
            xlim([0, 1])
            ylim([0, 1])
            set(gca, 'ydir', 'reverse')
            legend('Left Eye', 'Right Eye')
        end
        
        function val = isempty(obj)
            val = isempty(obj.X);
        end
%         
%         function val = length(obj)
%             warning('This has changed')
%             val = length(obj.X);
%         end
        
    end
    
    methods (Hidden)
        
%         function [vl, vr] = calcVelocity(obj)
%             ldx = [nan, diff(obj.LeftX)];
%             ldy = [nan, diff(obj.LeftY)];
%             rdx = [nan, diff(obj.RightX)];
%             rdy = [nan, diff(obj.RightY)];
%             vl = sqrt((ldx .^ 2) + (ldy .^ 2));
%             vr = sqrt((rdx .^ 2) + (rdy .^ 2));
%         end
%         
%         function acc = calcAcceleration(obj)
%             [vl, vr] = obj.calcVelocity;
%             
%         end
%         
%         function dir = calcDirection(obj)
%         end
%         
%         function [ld, rd] = calcDistances(obj)
%             ldx = [nan, diff(obj.LeftX)];
%             ldy = [nan, diff(obj.LeftY)];
%             rdx = [nan, diff(obj.RightX)];
%             rdy = [nan, diff(obj.RightY)];
%             ld = sqrt((ldx .^ 2) + (ldy .^ 2));
%             rd = sqrt((rdx .^ 2) + (rdy .^ 2));
%         end
            
        
    end
%     
%     methods (Hidden, Access = protected)
%         
% 
%         
%     end

end
    
    