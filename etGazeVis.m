classdef etGazeVis < handle
    
    properties
        Gaze etGazeData
        HiddenEvents = {...
            '*FRAME_CAL*',...
            };
        Colours 
%         LeftEyeColour = [066, 133, 244]
%         RightEyeColour = [052, 168, 083]
%         AverageColour = [213, 008, 000]    
%         LeftEyeVelocityColour = [220, 080, 080]
%         RightEyeVelocityColour = [220, 020, 020]
%         EventColour = [255, 245, 50]
        ShowMedianFilter = true;
    end
    
    properties (Dependent)
        ViewportEdge 
        ViewportWidth
        Cursor
    end 
    
    properties (Dependent, SetAccess = private)
        Duration
        DurationFormatted
    end    
    
    properties (Access = private)
        prViewportEdge = 0
        prViewportWidth = 10
        prCursor
        prVisibleGaze
        prVisibleTime
        prEventLabelExtent 
        uiFig
        uiX
        uiY
        uiToolbar
        uiBtnZoomOut
        uiBtnZoomIn
        uiBtnPageBack
        uiBtnPageForward
        uiCursorX
        uiCursorY
        filtMedian
        
        uiCol
    end
    
    methods
        
        function obj = etGazeVis(varargin)
            
            % store gaze data, if supplied
            if nargin == 1 && isa(varargin{1}, 'etGazeData')
                obj.Gaze = varargin{1};
            end
            
            % scale colours to normalised range
            obj.Colours = teCollection;
            obj.Colours('LeftEye') =                (   [066, 133, 244]     ./ 255);
            obj.Colours('RightEye') =               (   [052, 168, 083]     ./ 255);
            obj.Colours('Average') =                (   [213, 008, 000]     ./ 255);
            obj.Colours('LeftVelocity') =           (   [220, 180, 180]     ./ 255);
            obj.Colours('RightVelocity') =          (   [220, 020, 020]     ./ 255);
            obj.Colours('LeftAcceleration') =       (   [230, 230, 080]     ./ 255);
            obj.Colours('RightAcceleration') =      (   [220, 200, 020]     ./ 255);
            obj.Colours('LeftDirection') =          'jet';
            obj.Colours('RightDirection') =         'jet';
            obj.Colours('Events') =                 (   [255, 245, 050]     ./ 255);
            
            % setup UI
            obj.calcAxisLimits
            obj.uiCreate
            obj.uiDraw
        end
        
        function delete(obj)
            if ~isempty(obj.uiFig) && isgraphics(obj.uiFig)
                close(obj.uiFig)
            end
        end
        
        function Refresh(obj)
            obj.uiDraw
        end

        % get / set
        function val = get.Duration(obj)
            val = obj.Gaze.Duration;
        end
        
        function val = get.ViewportEdge(obj)
            val = obj.prViewportEdge;
        end
        
        function set.ViewportEdge(obj, val)
            if val > obj.Duration
                val = obj.Duration - obj.ViewportWidth;
            elseif val < 0
                val = 0;
            end
            obj.prViewportEdge = val;
            obj.calcAxisLimits
            obj.uiDraw 
        end
        
        function val = get.ViewportWidth(obj)
            val = obj.prViewportWidth;
        end
        
        function set.ViewportWidth(obj, val)
            if val < .5
                val = 0.5;
            elseif val > 600
                val = 600;
            end
            obj.prViewportWidth = val;
            obj.calcAxisLimits
            obj.uiDraw
        end
        
        function val = get.DurationFormatted(obj)
            val = datestr(obj.Duration / 86400, 'HH:MM:SS.fff');
        end
        
        function val = get.Cursor(obj)
            val = obj.prCursor;
        end
        
        function set.Cursor(obj, val)
            if val < 0
                error('Cursor cannot be set to a negative value.')
            elseif val > obj.Duration 
                error('Cursor cannot be set to a value greater than the total duration (%.2f)',...
                    obj.Duration)
            end
            obj.prCursor = val;
            obj.updateCursor
        end
        
        function set.prVisibleGaze(obj, val)
            obj.prVisibleGaze = val;
            obj.updateTimeSeriesGazeTime
        end
        
        function set.prVisibleTime(obj, val)
            obj.prVisibleTime = val;
            obj.updateTimeSeriesGazeTime
        end
        
    end
    
    methods (Hidden)
        
        % UI
        function uiCreate(obj)
            
            monitorPos = get(groot, 'MonitorPositions');
            numMonitors = size(monitorPos, 1); 
            if numMonitors > 1
                % draw figure full screen on last monitor
                rect = monitorPos(numMonitors, :);
                windowStyle = 'normal';
            else
                % draw figure docked
                rect = [0, 0, 100, 100];
                windowStyle = 'docked';
            end
            
            obj.uiFig = figure(...
                'visible', 'off',...
                'MenuBar', 'none',...
                'units', 'pixels',...
                'position', rect,...
                'windowStyle', windowStyle,...
                'SizeChangedFcn', @obj.uiResize);
            whitebg(obj.uiFig, 'k')
            
            pos = obj.uiGetPositions;
            
            obj.uiCol = etTimeseriesCollection('Position', pos.uiAxes);
            obj.uiCol('x') = etTimeseries(obj.prVisibleGaze, {'LeftX', 'RightX'},...
                'Time', obj.prVisibleTime,...
                'Parent', obj.uiFig,...
                'YLabel', 'Screen Coords',...
                'Colour', {obj.Colours('LeftEye'), obj.Colours('RightEye')});
            obj.uiCol('y') = etTimeseries(obj.prVisibleGaze, {'LeftY', 'RightY'},...
                'Time', obj.prVisibleTime,...
                'Parent', obj.uiFig,...
                'YLabel', 'Screen Coords',...
                'Colour', {obj.Colours('LeftEye'), obj.Colours('RightEye')});
            
            obj.uiCol('vel') = etTimeseries(obj.prVisibleGaze, {'Calc.VelocityLeft', 'Calc.VelocityRight'},...
                'Time', obj.prVisibleTime,...
                'YDir', 'normal',...                
                'Parent', obj.uiFig,...
                'YLabel', 'Screen Coords',...
                'Colour', {obj.Colours('LeftVelocity'), obj.Colours('RightVelocity')});
            
            obj.uiCol('acc') = etTimeseries(obj.prVisibleGaze, {'Calc.AccelerationLeft', 'Calc.AccelerationRight'},...
                'Time', obj.prVisibleTime,...
                'YDir', 'normal',...
                'Parent', obj.uiFig,...
                'YLabel', 'Screen Coords',...
                'Colour', {obj.Colours('LeftAcceleration'), obj.Colours('RightAcceleration')});  
            
            obj.uiCol('dir') = etTimeseries(obj.prVisibleGaze, {'Calc.DirectionLeft', 'Calc.DirectionRight'},...
                'Time', obj.prVisibleTime,...
                'YDir', 'normal',...
                'YLim', [0, 360],...
                'Parent', obj.uiFig,...
                'YLabel', '°',...
                'Colour', {obj.Colours('LeftDirection'), obj.Colours('RightDirection')});               
            
            %             obj.uiCreateTimeseries
            
%             obj.uiX = axes(...
%                 'units', 'pixels',...
%                 'position', pos.uiX,...
%                 'color', [0.1, 0.1, 0.1],...
%                 'pickableparts', 'visible',...
%                 'hittest', 'on',...
%                 'buttondownfcn', @obj.uiXY_Click);
%             
%             obj.uiY = axes(...
%                 'units', 'pixels',...
%                 'position', pos.uiY,...
%                 'color', [0.1, 0.1, 0.1],...
%                 'pickableparts', 'visible',...
%                 'hittest', 'on',...
%                 'buttondownfcn', @obj.uiXY_Click);
            
            obj.uiToolbar = uipanel('units', 'pixels', 'position', pos.uiToolbar);
            
            obj.uiBtnZoomOut = uicontrol(...
                'parent', obj.uiToolbar,...
                'units', 'pixels',...
                'position', pos.uiBtnZoomOut,...
                'style', 'pushbutton',...
                'string', '-',...
                'callback', @obj.uiBtnZoomOut_Click);
            
            obj.uiBtnZoomIn = uicontrol(...
                'parent', obj.uiToolbar,...            
                'units', 'pixels',...
                'position', pos.uiBtnZoomIn,...
                'style', 'pushbutton',...
                'string', '+',...
                'callback', @obj.uiBtnZoomIn_Click);
            
            obj.uiBtnPageBack = uicontrol(...
                'parent', obj.uiToolbar,...            
                'units', 'pixels',...
                'position', pos.uiBtnPageBack,...
                'style', 'pushbutton',...
                'string', '<',...
                'callback', @obj.uiBtnPageBack_Click);          
            
            obj.uiBtnPageForward = uicontrol(...
                'parent', obj.uiToolbar,...            
                'units', 'pixels',...
                'position', pos.uiBtnPageForward,...
                'style', 'pushbutton',...
                'string', '>',...
                'callback', @obj.uiBtnPageForward_Click);              
            
            obj.uiFig.Visible = 'on';
                        
        end
        
        function uiCreateTimeseries(obj)

        end
        
        function uiResize(obj, ~, ~)
            
            if ~isequal(obj.uiFig.Visible, 'on')
                return
            end
            
            pos = obj.uiGetPositions;
            
            obj.uiCol.Position = pos.uiAxes;
            obj.uiX.Position = pos.uiX;
            obj.uiY.Position = pos.uiY;
            obj.uiToolbar.Position = pos.uiToolbar;
            obj.uiBtnZoomOut = pos.uiBtnZoomOut;
            obj.uiBtnZoomIn = pos.uiBtnZoomIn;
            obj.uiBtnPageBack = pos.uiBtnPageBack;
            obj.uiBtnPageForward = pos.uiBtnPageForward;
                        
        end
        
        function uiDraw(obj)
            
            if isempty(obj.uiCol)
                return
            end
            
            for i = 1:obj.uiCol.Count
                item = obj.uiCol.Items{i};
                item.Gaze = obj.prVisibleGaze;
                item.Time = obj.prVisibleTime;
                item.Draw;
            end
            
%             cla(obj.uiX)
%             cla(obj.uiY)
%             hold(obj.uiX, 'on');
%             hold(obj.uiY, 'on');
%             
%             obj.plotMissingData
%             obj.plotEvents
%             obj.plotGaze
%             obj.drawCursor

        end
        
        function pos = uiGetPositions(obj)
            
            w = obj.uiFig.Position(3);
            h = obj.uiFig.Position(4);
            
            h_toolbar = 40;
            
            w_button = h_toolbar - 4;
            h_button = w_button;
            yoff_button = 2;
            
            h_axis = (h - h_button) / 2;
            yoff_axis = 20;
            
            pos.uiAxes = [0, h_toolbar, w, h - h_toolbar];
            
            pos.uiX = [0, h_toolbar + h_axis, w, h_axis - yoff_axis];
            pos.uiY = [0, h_toolbar + yoff_axis, w, h_axis - yoff_axis];
            pos.uiToolbar = [0, 0, w, h_toolbar];
            
            pos.uiBtnZoomOut = [0 * w_button, yoff_button, w_button, h_toolbar - (yoff_button * 2)];
            pos.uiBtnZoomIn = [1 * w_button, yoff_button, w_button, h_toolbar - (yoff_button * 2)];
            pos.uiBtnPageBack = [2 * w_button, yoff_button, w_button, h_toolbar - (yoff_button * 2)];
            pos.uiBtnPageForward = [3 * w_button, yoff_button, w_button, h_toolbar - (yoff_button * 2)];
            
        end
        
        % plotting
        function plotMissingData(obj)
            
%             % plot missing data
%             ctt = contig2time(findcontig2(obj.prVisibleGaze.Missing, true), obj.prVisibleTime);
%             for i = 1:size(ctt, 1)
%                 x1 = ctt(i, 1);
%                 x2 = ctt(i, 3);
%                 y1 = 0.00;
%                 y2 = 1.00;
%                 
%                 % X missing BG red
%                 rectangle(obj.uiX,...
%                     'position', [x1, y1, x2, y2],...
%                     'FaceColor', [1, 0, 0, .025],...
%                     'EdgeColor', 'none',...
%                     'PickableParts', 'all',...
%                     'HitTest', 'off');
% 
%                 % Y missing BG red
%                 rectangle(obj.uiY,...
%                     'position', [x1, y1, x2, y2],...
%                     'FaceColor', [1, 0, 0, .025],...
%                     'EdgeColor', 'none',...
%                     'PickableParts', 'all',...
%                     'HitTest', 'off');
%                 
%                 % Y missing lower red highlight
%                 missingRect = rectangle(obj.uiY,...
%                     'position', [x1, 0.98, x2, 0.02],...
%                     'FaceColor', [1, 0, 0, .75],...
%                     'EdgeColor', 'none');
%                 
%             end
            
        end
        
        function plotEvents(obj)
            
%             % plot events
%             tx = [];
%             tb = [];
%             if ~isempty(obj.prVisibleGaze.Events)
%                 ty = 0;
%                 for i = 1:size(obj.prVisibleGaze.Events, 1)
%                     
%                     x = obj.prVisibleGaze.Events.Time(i);
%                     lab = obj.prVisibleGaze.Events.Label{i};
%                     
%                     line(obj.uiX, [x, x], [0, 1],...
%                         'Color', [obj.EventColour, 0.5],...
%                         'linewidth', 3,...
%                         'hittest', 'off')
%                     
%                     line(obj.uiY, [x, x], [0, 1],...
%                         'Color', [obj.EventColour, 0.5],...
%                         'linewidth', 3,...
%                         'hittest', 'off')
%                     
%                     if ~isempty(tx)
%                         tb = tx.Extent;
%                         moveDown = x < tb(1) + tb(3);
%                     else
%                         moveDown = false;
%                     end
%                     
%                     if moveDown
%                         ty = ty + tb(4);
%                     else
%                         ty = 0;
%                     end
%                     
%                     tx = text(obj.uiX, x, ty, lab,...
%                         'Interpreter', 'none',...
%                         'Color', 'w',...
%                         'PickableParts', 'all',...
%                         'HitTest', 'off');
%                     
%                 end
%             end
            
        end
        
        function plotGaze(obj)
            
            markerAlpha = 0.5;
            
            % plot x
            sc_lx = scatter(obj.prVisibleTime, obj.prVisibleGaze.LeftX,...
                'parent', obj.uiX,...
                'pickableparts', 'all',...
                'hittest', 'off');
            sc_lx.MarkerEdgeColor = obj.Colours('LeftEye');
            sc_lx.MarkerFaceColor = sc_lx.MarkerEdgeColor;
            sc_lx.MarkerFaceAlpha = markerAlpha;     
            
            sc_rx = scatter(obj.prVisibleTime, obj.prVisibleGaze.RightX,...
                'parent', obj.uiX,...
                'pickableparts', 'all',...
                'hittest', 'off');
            sc_rx.MarkerEdgeColor = obj.Colours('RightEye');
            sc_rx.MarkerFaceColor = sc_rx.MarkerEdgeColor;
            sc_rx.MarkerFaceAlpha = markerAlpha; 
            
            obj.uiX.XGrid = 'on';
            obj.uiX.XMinorGrid = 'on';
            xlim(obj.uiX, [obj.ViewportEdge, obj.ViewportEdge + obj.ViewportWidth]);
            ylim(obj.uiX, [0, 1])
            obj.uiX.YDir = 'reverse';
            
            if obj.ShowMedianFilter 
                sc_lx_mf = plot(obj.prVisibleTime, obj.filtMedian.LeftX,...
                    'parent', obj.uiX,...
                    'color', 'w',...
                    'linewidth', 2,...
                    'pickableparts', 'all',...
                    'hittest', 'off');
                sc_rx_mf = plot(obj.prVisibleTime, obj.filtMedian.RightX,...
                    'parent', obj.uiX,...               
                    'color', 'w',...
                    'linewidth', 2,...
                    'pickableparts', 'all',...
                    'hittest', 'off');          
            end
                        
            % plot y
            sc_ly = scatter(obj.prVisibleTime, obj.prVisibleGaze.LeftY,...
                'parent', obj.uiY,...
                'pickableparts', 'all',...
                'hittest', 'off');
            sc_ly.MarkerEdgeColor = obj.Colours('LeftEye');
            sc_ly.MarkerFaceColor = sc_ly.MarkerEdgeColor;
            sc_ly.MarkerFaceAlpha = markerAlpha;           
            
            sc_ry = scatter(obj.prVisibleTime, obj.prVisibleGaze.RightY,...
                'parent', obj.uiY,...
                'pickableparts', 'all',...
                'hittest', 'off');
            sc_ry.MarkerEdgeColor = obj.Colours('RightEye');
            sc_ry.MarkerFaceColor = sc_ry.MarkerEdgeColor;
            sc_ry.MarkerFaceAlpha = markerAlpha;              
            
            obj.uiY.XGrid = 'on';
            obj.uiY.XMinorGrid = 'on';      
            xlim(obj.uiY, [obj.ViewportEdge, obj.ViewportEdge + obj.ViewportWidth]);
            ylim(obj.uiY, [0, 1])
            obj.uiY.YDir = 'reverse';    
            
            obj.uiY.XTickLabels = cellfun(@(x)...
                datestr(str2double(x) / 86400, 'HH:MM:SS.fff'),...
                obj.uiY.XTickLabels, 'UniformOutput', false);
            
            if obj.ShowMedianFilter
                sc_ly_mf = plot(obj.prVisibleTime, obj.filtMedian.LeftY,...
                    'parent', obj.uiY,...
                    'color', 'w',...
                    'linewidth', 2,...
                    'pickableparts', 'all',...
                    'hittest', 'off');   
                sc_ry_mf = plot(obj.prVisibleTime, obj.filtMedian.RightY,...
                    'parent', obj.uiY,...               
                    'color', 'w',...
                    'linewidth', 2,...
                    'pickableparts', 'all',...
                    'hittest', 'off');       
            end
            
            legend([sc_lx, sc_rx], 'Left', 'Right')
            
        end
        
        function drawCursor(obj)
            
%             if isempty(obj.prCursor)
%                 return
%             end
%             
%             x = obj.Cursor;
%             
%             if isempty(obj.uiCursorX) || ~isgraphics(obj.uiCursorX)
%                 obj.uiCursorX = line(obj.uiX, [x, x], [0, 1],...
%                     'Color', 'm',...
%                     'linewidth', 2);
%             end
%                 
%             if isempty(obj.uiCursorY) || ~isgraphics(obj.uiCursorY)
%                 obj.uiCursorY = line(obj.uiY, [x, x], [0, 1],...
%                     'Color', 'm',...
%                     'linewidth', 2);   
%             end
            
        end
        
        function updateCursor(obj)
            
%             if isempty(obj.prCursor)
%                 return
%             end
%             
%             if isempty(obj.uiCursorX) || isempty(obj.uiCursorY)
%                 obj.drawCursor
%                 return
%             else
%                 obj.uiCursorX.XData = repmat(obj.Cursor, 1, 2);
%                 obj.uiCursorY.XData = repmat(obj.Cursor, 1, 2);
%             end
            
        end
        
        % callbacks
        function uiBtnZoomOut_Click(obj, ~, ~)
            obj.ViewportWidth = obj.ViewportWidth * 2;
        end
        
        function uiBtnZoomIn_Click(obj, ~, ~)
            obj.ViewportWidth = obj.ViewportWidth / 2;
        end
        
        function uiBtnPageBack_Click(obj, ~, ~)
            obj.ViewportEdge = obj.ViewportEdge - obj.ViewportWidth;      
            obj.uiDraw
        end
        
        function uiBtnPageForward_Click(obj, ~, ~)
            newEdge = obj.ViewportEdge + obj.ViewportWidth;
            if newEdge > obj.Gaze.Duration - obj.ViewportWidth
                newEdge = obj.Gaze.Duration - obj.ViewportWidth;
            end
            obj.ViewportEdge = newEdge;
            obj.uiDraw
        end
        
        function uiXY_Click(obj, ~, event)
            
            obj.Cursor = obj.ViewportEdge + event.IntersectionPoint(1);
            
        end

        % utils
        function calcAxisLimits(obj)
        % called when either the ViewportEdge position or the visible duration 
        % (zoom) changes. Selects the relevant portion of data from the
        % gaze object and stores it in prVisibleGaze property. From here,
        % the actual axes limits will be automatically set when the gaze
        % data is drawn
        
            % we will segment between times t1 and t2. t1 is the lefthand
            % edge of the viewport
            t1 = obj.prViewportEdge;
            
            % t2 is the righthand edge of the viewport. Now segment the
            % gaze data
            t2 = t1 + obj.prViewportWidth;
            
            % the currently set viewport width cannot exceed the duration
            % of gaze data, correct it if it does. Basically means you
            % cannot scroll/zoom past the edge of the data
            if t1 + t2 > obj.Gaze.Duration
                t2 = obj.Gaze.Duration;
                % update viewport width
                obj.prViewportWidth = t2 - t1;
            end
            
            obj.prVisibleGaze = obj.Gaze.SegmentByTime(t1, t2);
            
            % the timestamps in the segmented gaze data will be zeroed. So
            % add the lefthand viewport edge (aka t1) to each gaze sample's
            % timestamp to correct it. 
            obj.prVisibleTime = obj.prViewportEdge + obj.prVisibleGaze.Time;
            
            % remove any events in the .HiddenEvents property
            obj.filterOutHiddenEvents
            
            % (optionally) apply a filter to the gaze data 
            obj.filterGaze
            
        end
        
        function updateTimeSeriesGazeTime(obj)
            if isempty(obj.uiCol), return, end
            for i = 1:obj.uiCol.Count
                item = obj.uiCol.Items{i};
                item.Gaze = obj.prVisibleGaze;
                item.Time = obj.prVisibleTime;
            end
        end
        
        function filterOutHiddenEvents(obj)
            
            % if no events, nothing to filter
            if isempty(obj.prVisibleGaze.Events), return, end            
            
            ev = obj.prVisibleGaze.Events(:, 2);
            ev_char = cellfun(@cell2char, ev.Label, 'uniform', false);
            
            numFilters = length(obj.HiddenEvents);
            idx = false(size(ev_char));
            for f = 1:numFilters
                mask = regexptranslate('wildcard', obj.HiddenEvents{f});
                res = regexp(ev_char, mask);
                idx = idx | cellfun(@(x) ~isempty(x) && x == 1, res);                
            end
            obj.prVisibleGaze.Events(idx, :) = [];
            
        end                
        
        function filterGaze(obj)
            
            if obj.ShowMedianFilter
                obj.applyMedianFilter
            else
                obj.clearMedianFilter
            end
            
        end
        
        function calcEventLabelExtent(obj)
            tx = text(obj.uiX, 0, 0, 'TMP');
            obj.prEventLabelExtent = tx.Extent;
            delete(tx)
        end
        
        % filtering
        function applyMedianFilter(obj)
            % todo - should this be done on the main gaze object, just
            % once?
            % todo - handle weird sample rate issue when lots of missing
            % gaze
            span_sec = 0.1;
            span_samp = round(span_sec * obj.prVisibleGaze.MeanSampleRate);
            if isnan(span_samp), return, end
            raw = [...
                obj.prVisibleGaze.LeftX,...
                obj.prVisibleGaze.LeftY,...
                obj.prVisibleGaze.RightX,...
                obj.prVisibleGaze.RightY,...
                ];
            filt = medfilt1(raw, span_samp);
            obj.filtMedian.LeftX = filt(:, 1);
            obj.filtMedian.LeftY = filt(:, 2);
            obj.filtMedian.RightX = filt(:, 3);
            obj.filtMedian.RightY = filt(:, 4);
        end
        
        function clearMedianFilter(obj)
            obj.filtMedian = [];
        end
        
    end
        
end
   