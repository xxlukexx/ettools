function [area, tab] = etCalculateAOIArea(path_aoi, def, maxDur)

    clear port aoi tl
    sca

    port = teViewport;
%     load(path_def, 'def')
    
    aoi = vpaDynamicAOI;
    aoi.AOIDefinition = def;
    aoi.Path = path_aoi;
    port.Viewpane('aoi') = aoi;
    
    if ~exist('maxDur', 'var') || isempty(maxDur)
        maxDur = aoi.Duration;
    end
    tl = vpaTimeline;
    tl.Duration = maxDur;
    tl.ZPosition = inf;
    port.Viewpane('tl') = tl;
    aoi.Timeline = tl;

    [area, tab] = aoi.CalculateAOIArea(maxDur);
    tab = vertcat(tab{~cellfun(@isempty, tab)});

    [~, fil] = fileparts(path_aoi);
    figure('name', fil)
    for i = 1:size(area, 2)
        histogram(area(:, i), 'BinWidth', .02, 'DisplayStyle', 'stairs', 'LineWidth', 3)
        hold on
    end
    legend(def(:, 1))    
    
    clear port aoi tl
    sca
    
end