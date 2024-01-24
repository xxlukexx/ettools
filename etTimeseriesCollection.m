classdef etTimeseriesCollection < teCollection
    
    properties
        Position
        PropHeight
    end
    
    properties (Dependent)
        Order
    end
    
    properties (Access = private)
        prOrder
    end
    
    methods
        
        function obj = etTimeseriesCollection(varargin)
            obj = obj@teCollection('etTimeseries');
            obj.Position = varargin{2};
            % todo - input parser
        end
        
        function AddItem(obj, varargin)
            AddItem@teCollection(obj, varargin{:})
            obj.UpdateHeights
        end
        
        function UpdateHeights(obj)
            h = obj.Position(4);
            w = obj.Position(3);
            ydiv = h * 0.01;
            obj.PropHeight = repmat(1 / obj.Count, obj.Count, 1);
%             ph = obj.Items{1}.Parent.Position(4);
%             pw = obj.Items{1}.Parent.Position(3);
            y1 = obj.Position(2);
            items = flipud(obj.Items);
            for i = obj.Count:-1:1
                items{i}.Position = [35, y1 + ydiv, w - 40, (h * obj.PropHeight(i)) - (ydiv * 2)];
                y1 = y1 + (h * obj.PropHeight(i));
                items{i}.IsTop = false;
                items{i}.IsBottom = false;
                if i == obj.Count
                    items{i}.IsBottom = true;
                elseif i == 1
                    items{i}.IsTop = true;
                end
            end
            obj.Items = flipud(items);
        end
        
        function UpdateOrder(obj)
        end
        
        % get/set
        function val = get.Order(obj)
            val = obj.prOrder;
        end
        
        function set.Order(obj, val)
            % todo - check bounds, uniqueness etc
            obj.prOrder = val;
        end
        
        function set.PropHeight(obj, val)
            if isempty(val), return, end
            if ~isvector(val) || ~isnumeric(val) || length(val) ~= obj.Count
                error('PropHeight must be a numeric vector, with an element for each element in the collection.')
            end
            if ~sum(val)
                val(end) = 1 - val(1:end - 1);
                warning('Adjusted final height to make proportions sum to 1.')
            end
            obj.PropHeight = val;
        end
        
        function set.Position(obj, val)
            obj.Position = val;
            obj.UpdateHeights
        end
        
    end
    
end
                
            
        

    