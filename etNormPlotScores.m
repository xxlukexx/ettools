function etNormPlotScores(scores, gaze)

    if ~iscell(scores)
        scores = {scores};
    end
    
    if ~iscell(gaze)
        gaze = {gaze};
    end
    
    figure
    
    numGroups = length(scores);
    for g = 1:numGroups
        
        mu = nanmean(scores{g}, 2)';
        sd = nanstd(scores{g}, [], 2)';

        mu = smooth(mu, 100);
        sd = smooth(sd, 100);

        plot(gaze{g}.Time, mu, 'LineWidth', 3)
        hold on

        set(gca, 'ColorOrderIndex', g)
        plot(gaze{g}.Time, mu + sd, 'LineStyle', '--')
        set(gca, 'ColorOrderIndex', g)
        plot(gaze{g}.Time, mu - sd, 'LineStyle', '--')

        set(gca, 'ColorOrderIndex', g)
        plot(gaze{g}.Time, mu + (2 * sd), 'LineStyle', '-.')
        set(gca, 'ColorOrderIndex', g)
        plot(gaze{g}.Time, mu - (2 * sd), 'LineStyle', '-.')
        
    end
    
end
    


