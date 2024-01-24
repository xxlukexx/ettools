classdef etTimeseries < dynamicprops
    
    properties
        Gaze
        Time
        Property
        Parent
        Colour
        Visible 
        ShowMissingData = true
        ShowEvents = true
        IsTop = false
        IsBottom = false
        H
        EventColour = [1.00, 0.96, 0.20]
        YDir = 'reverse'
        YLim = [0, 1]
        YLabel = []
    end
    
    properties (Dependent)
        Position
        Cursor
    end
    
    properties (Dependent, SetAccess = private)
        NumProperties = 0
    end
    
    properties (Access = private)
        prFigureOwner = false
        prCursor
        % UI
        uiCursorX
        uiCursorY
        uiUpdating = false
    end
    
    methods
        
        function obj = etTimeseries(gaze, prop, varargin)
            
            parser         =   inputParser;
            isGazeObj      =   @(x) isa(x, 'etGazeData');
            isTimeVec      =   @(x) isvector(x) && isnumeric(x);
            
            addRequired(    parser, 'gaze',                     @(x) isa(x, 'etGazeData'))
            addRequired(    parser, 'property',                 @(x) ischar(x) || iscellstr(x))
            
            addParameter(   parser, 'time',                     [],     @(x) isvector(x) && isnumeric(x))
            addParameter(   parser, 'parent',                   [],     @isgraphics )
            addParameter(   parser, 'position',                 [],     @(x) isvector(x) && length(x) == 4)
            addParameter(   parser, 'visible',                  true,   @islogical)
            addParameter(   parser, 'colour',                   lines(1),@(x) iscell(x) || isvector(x))
            addParameter(   parser, 'ydir',                     'reverse', @ischar)
            addParameter(   parser, 'ylim',                     [0, 1])
            addParameter(   parser, 'ylabel',                   [])
            
            parse(          parser, gaze, prop, varargin{:});
            obj.Gaze        =   parser.Results.gaze;
            obj.Property    =   parser.Results.property;
            obj.Time        =   parser.Results.time;
            obj.Position    =   parser.Results.position;
            obj.Parent      =   parser.Results.parent;
            obj.Visible     =   parser.Results.visible;
            obj.Colour      =   parser.Results.colour;
            obj.YDir        =   parser.Results.ydir;
            obj.YLim        =   parser.Results.ylim;
            obj.YLabel      =   parser.Results.ylabel;
            
            allProps = properties(obj.Gaze);
            calcProps = fieldnames(obj.Gaze.Calc);
            allProps = [allProps; calcProps];
            flatProps = strrep(obj.Property, 'Calc.', '');
            if ~all(cellfun(@(x) ismember(x, allProps), flatProps))
                error('The .Property function (%s) is not a valid property of an etGazeData object.',...
                    obj.Property)
            end
            
            if isempty(obj.Parent)
                obj.Parent = figure(...
                    'units',            'pixels',...
                    'SizeChangedFcn',   @obj.uiFigureResize);
                obj.prFigureOwner = true;
            end
            
            obj.uiCreate
            
        end
        
        function delete(obj)
            if obj.prFigureOwner
                delete(obj.Parent)
            end
        end
        
        function Draw(obj)
            
            % if the UI is already updating, do not attempt another update
            if obj.uiUpdating, return, end
            
            % set the updating flag to prevent further updates until we are
            % done here
            obj.uiUpdating = true;
            
            % clear the axes and draw all components
            cla(obj.H)
            hold(obj.H, 'on')
            obj.plotMissingData
            obj.plotEvents
            obj.plotGaze
            obj.drawCursor
            
            % unset the updating flag to allow future updates
            obj.uiUpdating = false;
            
        end
        
        function plotMissingData(obj)
            
            if ~obj.ShowMissingData, return, end
            
            if ~isempty(obj.Time)
                t = obj.Time;
            else
                t = obj.Gaze.Time;
            end
            
            % plot missing data
            ctt = contig2time(findcontig2(obj.Gaze.Missing, true), t);
            for i = 1:size(ctt, 1)
                x1 = ctt(i, 1);
                x2 = ctt(i, 3);
                y1 = 0.00;
                y2 = 1.00;

                % missing BG red
                rectangle(obj.H,...
                    'position', [x1, y1, x2, y2],...
                    'FaceColor', [1, 0, 0, .05],...
                    'EdgeColor', 'none',...
                    'PickableParts', 'all',...
                    'HitTest', 'off');
                
                if obj.IsBottom
                    % missing lower red highlight
                    rectangle(obj.H,...
                        'position', [x1, 0.98, x2, 0.02],...
                        'FaceColor', [1, 0, 0, .75],...
                        'EdgeColor', 'none')
                end
                
            end
            
        end
        
        function plotEvents(obj)
            
            % plot events
            tx = [];
            tb = [];
            if ~isempty(obj.Gaze.Events)
                ty = 0;
                for i = 1:size(obj.Gaze.Events, 1)
                    
                    x = obj.Gaze.Events.Time(i);
                    lab = obj.Gaze.Events.Label{i};
                    
                    line(obj.H, [x, x], [0, 1],...
                        'Color', [obj.EventColour, 0.5],...
                        'linewidth', 3,...
                        'hittest', 'off')
                    
                    if ~isempty(tx)
                        tb = tx.Extent;
                        moveDown = x < tb(1) + tb(3);
                    else
                        moveDown = false;
                    end
                    
                    if moveDown
                        ty = ty + tb(4);
                    else
                        ty = 0;
                    end
                    
                    % labels
                    if obj.IsTop
                        tx = text(obj.H, x, ty, lab,...
                            'Interpreter', 'none',...
                            'Color', 'w',...
                            'PickableParts', 'all',...
                            'HitTest', 'off');
                    end
                    
                end
            end
            
        end
        
        function plotGaze(obj)
            
            markerAlpha = 0.8;
            
            if ~isempty(obj.Time)
                x = obj.Time;
            else
                x = obj.Gaze.Time;
            end
%             x = obj.Gaze.Time;
            
            sc = {};
            for p = 1:obj.NumProperties
                
%                 y = obj.Gaze.(obj.Property{p});
                [y, col] = obj.calcY(obj.Property{p});
                
                sc{p} = scatter(x, y, [], col,...
                    'parent', obj.H,...
                    'pickableparts', 'all',...
                    'hittest', 'off');

                sc{p}.MarkerFaceColor = sc{p}.MarkerEdgeColor;
                sc{p}.MarkerFaceAlpha = markerAlpha;    
            
            end
            
            obj.H.XGrid = 'on';
            obj.H.XMinorGrid = 'on';
            xlim(obj.H, [x(1), x(end)]);
            ylim(obj.H, obj.YLim)
            obj.H.YDir = obj.YDir;
            
            legend([sc{:}], obj.Property)
            ylabel(obj.H, obj.YLabel)
            
        end
        
        function drawCursor(obj)
            
            if isempty(obj.prCursor)
                return
            end
            
            x = obj.Cursor;
            
            if isempty(obj.uiCursorX) || ~isgraphics(obj.uiCursorX)
                obj.uiCursorX = line(obj.uiX, [x, x], [0, 1],...
                    'Color', 'm',...
                    'linewidth', 2);
            end
                
            if isempty(obj.uiCursorY) || ~isgraphics(obj.uiCursorY)
                obj.uiCursorY = line(obj.uiY, [x, x], [0, 1],...
                    'Color', 'm',...
                    'linewidth', 2);   
            end
            
        end
        
        function updateCursor(obj)
            
            if isempty(obj.prCursor)
                return
            end
            
            if isempty(obj.uiCursorX) || isempty(obj.uiCursorY)
                obj.drawCursor
                return
            else
                obj.uiCursorX.XData = repmat(obj.Cursor, 1, 2);
                obj.uiCursorY.XData = repmat(obj.Cursor, 1, 2);
            end
            
        end
        
        % get/set
        function set.Position(obj, val)
            obj.H.Position = val;
        end
        
        function val = get.Position(obj)
            val = obj.H.Position;
        end
        
        function set.Property(obj, val)
            % todo - check this is a valid cellstr and each element is a
            % valid prop
            if ~iscell(val)
                val = {val};
            end
            obj.Property = val;
        end
        
        function val = get.NumProperties(obj)
            val = length(obj.Property);
        end
        
        function val = get.Cursor(obj)
            val = obj.prCursor;
        end
        
        function set.Cursor(obj, val)
            obj.prCursor = val;
        end
        
%         function set.IsTop(obj, val)
%             obj.IsTop = val;
%         end
        

        
        
        
        
        
%         function obj = etTimeseries(gaze, time, prop, fcn, colour)
%             isfun           =   @(f) isa(f, 'function_handle') || exist(f, 'file') == 2;
    end
    
    methods (Hidden)
        
        function uiCreate(obj)
            
            if isempty(obj.Position)
                obj.Position = [50, 0, obj.Parent.Position(3) - 50, obj.Parent.Position(4)];
            end
            
            obj.H = axes(...
                'units', 'pixels',...
                'outerposition', obj.Position,...
                'color', [0.1, 0.1, 0.1],...
                'pickableparts', 'visible',...
                'hittest', 'on',...
                'buttondownfcn', @obj.Click);   
            
        end
        
        function Click(~, ~, ~)
        end
        
        function uiFigureResize(obj, src, ~)
            if isempty(obj.Parent), return, end
            obj.Position = [0, 0, obj.Parent.Position(3:4)];
        end
        
        function [y, col] = calcY(obj, prop)
            if ~exist('prop', 'var') || isempty(prop)
                error('Must pass a property to calculate y for.')
            end
            if contains(prop, 'Calc.')
                prop = strrep(prop, 'Calc.', '');
                if ~hasField(obj.Gaze.Calc, prop)
                    error('Field not found in .Calc structure: %s', prop)
                end
                y = obj.Gaze.Calc.(prop);
            else
                if ~isprop(obj.Gaze, prop)
                    error('Property not found in Gaze object: %s', prop)
                end
                y = obj.Gaze.(prop);
            end
            allProps = strrep(obj.Property, 'Calc.', '');
            idx = find(cellfun(@(x) strcmpi(x, prop), allProps));
            if ischar(obj.Colour{idx})
                colMap = eval(sprintf('%s(length(y))', obj.Colour{idx}));
                rnk = tiedrank(y);
                col = nan(length(y), 3);
                missing = isnan(y);
                col(~missing, :) = colMap(rnk(~missing), :);
            else
                col = obj.Colour{idx};
            end
        end

    
    end
    
end
    
    