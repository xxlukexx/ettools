function [tab_stats, fig_hist, fig_boxplot, res_t] =...
    etPlotAOIFactorial(tab, measure, varargin)

    set(groot, 'defaultAxesTickLabelInterpreter','none')
    set(groot, 'defaultLegendInterpreter','none');
    set(groot, 'defaulttextInterpreter', 'none');

    parser      =   inputParser;
    checkField  =   @(x) any(strcmpi(tab.Properties.VariableNames, x)   );
    addRequired(    parser, 'tab',                              @istable        )
    addRequired(    parser, 'measure',                          checkField      )
    addParameter(   parser, 'compare',              []                          )
    addParameter(   parser, 'rows',                 [],         checkField      )
    addParameter(   parser, 'cols',                 [],         checkField      )
    addParameter(   parser, 'stack',                [],         checkField      )
    addParameter(   parser, 'name',                 'Average',  @ischar         )
    addParameter(   parser, 'colMap',               @lines                      )
    addParameter(   parser, 'linewidth',            1.5,        @isnumeric      )
    addParameter(   parser, 'fontsize',             10,         @isnumeric      )
    addParameter(   parser, 'plotBoxplot',          true,       @islogical      )
    addParameter(   parser, 'plotHist',             false,      @islogical      )
    addParameter(   parser, 'plotStacked',          false,      @islogical      )
    addParameter(   parser, 'makeTable',            false,      @islogical      )
    addParameter(   parser, 'legend',               true,       @islogical      )
    addParameter(   parser, 'ttest',                false,      @islogical      )
    addParameter(   parser, 'cohensd',              false,      @islogical      )
    addParameter(   parser, 'equalYLim',            false,      @islogical      )
    addParameter(   parser, 'yLabel',               [],         @ischar         )
    addParameter(   parser, 'xLabel',               [],         @ischar         )
    addParameter(   parser, 'smoothHist',           'none',     @ischar           )
    addParameter(   parser, 'colours',              [],         @isnumeric      )
    
    parse(          parser, tab, measure, varargin{:});
    tab         =   parser.Results.tab;
    measure     =   parser.Results.measure;
    compare     =   parser.Results.compare;
    rows        =   parser.Results.rows;
    cols        =   parser.Results.cols;
    stack       =   parser.Results.stack;
    name        =   parser.Results.name;
    colMap      =   parser.Results.colMap;
    linewidth   =   parser.Results.linewidth;
    fontsize    =   parser.Results.fontsize;
    plotBoxPlot =   parser.Results.plotBoxplot;
    plotHist    =   parser.Results.plotHist;
    plotStacked =   parser.Results.plotStacked;
    showLegend  =   parser.Results.legend;
    doTTest     =   parser.Results.ttest;
    doCohensD   =   parser.Results.cohensd;
    doTable     =   parser.Results.makeTable;
    equalYLim   =   parser.Results.equalYLim;
    yLabel      =   parser.Results.yLabel;
    xLabel      =   parser.Results.xLabel;
    smoothHist  =   parser.Results.smoothHist;
    colours     =   parser.Results.colours;
    
    coloursSpecified = ~isempty(colours);
    
    % y label by default is measure name, unless otherwise specified with
    % the yLabel option
    if isempty(yLabel)
        yLabel = measure;
    end
    if isempty(xLabel)
        xLabel = compare;
    end
    
    % histograms can be smoothed with 'spline', have a 'gaussian' fitted to
    % it
    if ~ismember(smoothHist, {'none', 'spline', 'gaussian'})
        error('The ''smoothHist'' parameter must be either ''none'', ''spline'' or ''gaussian''')
    end

    % axes line width is one thinner than line width 
    axLineWidth = linewidth - 1;
    if linewidth <= 2, axLineWidth = 2; end
    
    numData = size(tab, 1);
    
    % set up values of compare, stack, rows and cols variables
    [comp_u, comp_s, numComp] = prepareVariable(tab, compare);
    [row_u, row_s, numRow] = prepareVariable(tab, rows);
    [col_u, col_s, numCol] = prepareVariable(tab, cols);
    [stack_u, stack_s, numStack] = prepareVariable(tab, stack);
    
%     % cannot plots histograms if more than 2 comp
%     if numComp ~= 2 && plotB2BHist
%         warning('Back-to-back histograms can only be plotted when there are two comparisons.')
%         plotB2BHist = false;
%     end

    %% t-test
    
    p = [];
    mu = [];
    sd = [];
    df = [];
    tstat = [];
    if doTTest && numComp == 2
        
        % do t-test for each row/col
        p = nan(numRow, numCol);
        for c = 1:numCol
            for r = 1:numRow
                
                % filter table for this row, col
                idx_filt = row_s == r & col_s == c;
                tmp = tab(idx_filt, :);
                
                % make subs for comparison
                [comp_u, ~, comp_s] = unique(tmp.(compare));
                numComp = length(comp_u);
                
                x = tmp.(measure)(comp_s == 1);
                y = tmp.(measure)(comp_s == 2);
                
                [~, p(r, c), ~, stats_tmp] = ttest2(x, y);
                df(r, c) = stats_tmp.df;
                tstat(r, c) = stats_tmp.tstat;
%                 mu(r, c, 1) = nanmean(x);
%                 mu(r, c, 2) = nanmean(y);
%                 sd(r, c, 1) = nanmean(x);
%                 sd(r, c, 2) = nanmean(y);
                
                
            end
        end
        
        % convert p values to cell array of string
        t_cell = arrayfun(@(df, t, p) sprintf('t(%d)=%.2f, p=%.3f', df, t, p), df, tstat, p, 'uniform', false);
        
        % find any p values <.001 (which will display as p=0.000) and
        % replace with p<.001
        idx_low = p < .001;
%         t_cell(idx_low) = repmat({'p<.001'}, sum(idx_low(:), 1), 1);
        t_cell(idx_low) = strrep(t_cell(idx_low), 'p=0.000', 'p<.001');
        
        % format results in table
        tab_t = cell2table(t_cell, 'RowNames', row_u, 'VariableNames', col_u);
%         tab_mu = array2table(mu, 'RowNames', row_u, 'VariableNames', col_u);
%         tab_sd = array2table(sd, 'RowNames', row_u, 'VariableNames', col_u);
        fprintf('T-test on %s:\n\n', compare);
        disp(tab_t)        
        
%         fprintf('Means of %s:\n\n', compare);
%         disp(tab_mu)
%         fprintf('SD of %s:\n\n', compare);
%         disp(tab_sd)
        
    elseif doTTest && numComp ~= 2
        
        warning('Cannot perform t-test unless number of comparisons is equal to two.')           
    end

    %% cohen's d
    
    d = [];
    if doCohensD && numComp == 2
        
        % calculate d for each row/col
        d = nan(numRow, numCol);
        for c = 1:numCol
            for r = 1:numRow
                
                % filter table for this row, col
                idx_filt = row_s == r & col_s == c;
                tmp = tab(idx_filt, :);
                
                % make subs for comparison
                [comp_u, ~, comp_s] = unique(tmp.(compare));
                numComp = length(comp_u);

                % calculate mean and SEM for each comparison
                mu = accumarray(comp_s, tmp.(measure), [], @nanmean);
                sd = accumarray(comp_s, tmp.(measure), [], @nanstd); 
                
                % pooled SD
                x1 = tmp.(measure)(comp_s == 1);
                x2 = tmp.(measure)(comp_s == 2);
                [~, sdp] = sdpooled(x1, x2);
                
                % calculate N for each comp
%                 n = accumarray(comp_s, tmp.(measure), [], @(x) sum(~isnan(x)));
% 
%                 % pooled SD
% 
%                 % normalise N
%                 n = n / max(n);
% 
%                 % calculate pooled SD - weight SD by normalised N and take mean
% %                 sdp = mean(sd .* n);
% %                 sdp = mean(sd);
% %                 [sdp, ~] = sdpooled
%                 warning('Check calculation of pooled SD for cohen''s d!!')

                % calculate cohen's d
                d(r, c) = abs(diff(mu) / sdp);
                
            end
        end
        
        % format results in table
        tab_d = array2table(d, 'RowNames', row_u, 'VariableNames', col_u);
        fprintf('Cohen''s d:\n\n');
        disp(tab_d)
        
    elseif doCohensD && numComp ~= 2
        
            warning('Cannot calculate Cohen''s d unless number of comparisons is equal to two.')           
    end
   
    %% boxplots
    if plotBoxPlot
        
        fig_barhist = figure('name', name);
        set(fig_barhist, 'defaultaxesfontsize', fontsize);
        
        % set ylimits for all plots
        yMax = max(tab.(measure));
        yMin = min(tab.(measure));
        yl = [yMin, yMax];
        
        % loop through each comparison/row/col 
        spIdx = 1;
        for row = 1:numRow
            for col = 1:numCol
                
                % filter table for current row and col
                idx_filt = row_s == row & col_s == col;
                tmp = tab(idx_filt, :);
                
                % make subs for comparison
                [histcomp_u, ~, histcomp_s] = unique(tmp.(compare));
                
                % plot bar
                subplot(numRow, numCol, spIdx);
                
                notBoxPlot(tmp.(measure), histcomp_s, 'jitter', 0.5);
                set(gca, 'XTickLabel', histcomp_u);
                
                % make title
                if strcmpi(col_u{col}, 'NONE')
                    titleCol = '';
                else
                    titleCol = sprintf('%s: %s', cols, upper(col_u{col}));
                end
                if strcmpi(row_u{row}, 'NONE')
                    titleRow = '';
                else
                    titleRow = sprintf('%s: %s', rows, upper(row_u{row}));
                end
                titleStr = sprintf('%s %s', titleRow, titleCol);   
                
                if doCohensD
                    titleStr = sprintf('%s d=%.2f', d(row, col));
                end
                
                title(titleStr);
                if equalYLim, ylim(yl), end
                
                spIdx = spIdx + 1;
        
            end
            
        end

    end
    
    %% histograms
    if plotHist
        
        hist_smooth = 12;
        fig_barhist = figure('name', name);
        set(fig_barhist, 'defaultaxesfontsize', fontsize);
        
        % set ylimits for all plots
        yMax = max(tab.(measure));
        yMin = min(tab.(measure));
        yl = [inf, -inf];
                
        % loop through each comparison/row/col 
        spIdx = 1;
        sp = zeros(numRow * numCol, 1);
        yl_min = inf;
        yl_max = -inf;
        for row = 1:numRow
            for col = 1:numCol
                
                % filter table for current row and col
                idx_filt = row_s == row & col_s == col;
                tmp = tab(idx_filt,  :);
                
                % make subs for comparison
                [histcomp_u, ~, histcomp_s] = unique(tmp.(compare));
                histNumComp = length(histcomp_u);
                
                % calculate mean and SEM for each comparison
                bar_mu = accumarray(histcomp_s, tmp.(measure), [], @nanmean);
                if isempty(bar_mu), continue, end
                bar_sem = accumarray(histcomp_s, tmp.(measure), [], @nansem);
                
                % update ylimits
                if equalYLim
                    yl_new = [min(bar_mu - bar_sem), max(bar_mu + bar_sem)];
                    if yl_new(1) < yl(1), yl(1) = yl_new(1); end
                    if yl_new(2) > yl(2), yl(2) = yl_new(2); end
                end
                
                % plot bar
                sp(spIdx) = subplot(numRow, numCol * 2, spIdx);
                
                if ~coloursSpecified
                    colours = lines(histNumComp);
                else
                    set(gca, 'ColorOrder', colours)
                end
                
                % make title
                if strcmpi(col_u{col}, 'NONE')
                    titleCol = '';
                else
                    titleCol = sprintf('%s: %s', cols, upper(col_u{col}));
                end
                if strcmpi(row_u{row}, 'NONE')
                    titleRow = '';
                else
                    titleRow = sprintf('%s: %s', rows, upper(row_u{row}));
                end
                titleStr = sprintf('%s %s', titleRow, titleCol);
                
                % append t-test 
                if doTTest && ~isempty(p)
                    titleStr = sprintf('%s %s', titleStr, t_cell{row, col});
                end
                
                % append cohen's d
                if doCohensD && ~isempty(d)
                    titleStr = sprintf('%s d=%.2f', titleStr, d(row, col));
                end
                
                % plot bar
                [h_bar, h_errbar] = barwitherr(bar_sem * 2, bar_mu);
                h_bar.FaceColor = 'flat';
                for hc = 1:histNumComp
                    h_bar.CData(hc, :) = colours(hc, :);
                end
                h_bar.LineWidth = linewidth;
                set(h_errbar, 'LineWidth', linewidth);
                
                % settings
                xlabel(xLabel)
                ylabel(yLabel)
                set(gca, 'XTickLabel', comp_u)
                title(titleStr)
                ax = gca;
                ax.XRuler.Axle.LineWidth = axLineWidth;
                ax.YRuler.Axle.LineWidth = axLineWidth; 
                if ~any(isnan(bar_mu)) && ~any(isnan(bar_sem))
                    yl = [min(bar_mu - (bar_sem * 3)), max(bar_mu + (bar_sem *3))];
                end
                % the matchBarYLim flag means we want the same y axis
                % limits across all rows/cols. Otherwise we set the
                % ylim for the current row & col (i.e. current plot)
                % and leave it at that
                if equalYLim
                    % store the yl value if either edge exceeds the
                    % current range for all plots
                    if yl(1) < yl_min, yl_min = yl(1); end
                    if yl(2) > yl_max, yl_max = yl(2); end
                else
                    % set the ylim for this plot immediately, don't
                    % bother storing the value
                    ylim([yl(1) - abs(yl(1) * .2), yl(2) + abs(yl(2) * .2)])
                end
                
                spIdx = spIdx + 1;
 
                % plot hist
                subplot(numRow, numCol * 2, spIdx);
                hold on
                set(gca, 'ColorOrder', colours)
                
                for comp = 1:numComp
                    
                    % get data, make histogram
                    idx_comp = comp == histcomp_s;
                    idx_nan = isnan(tmp.(measure));
                    m = tmp.(measure)(idx_comp & ~idx_nan);
                    if isempty(m), continue, end
                    
                    hold on
                    
                    [vals, edges, ~] =...
                        histcounts(m, hist_smooth, 'Normalization', 'probability');                    
                    edges = edges(2:end) - (edges(2)-edges(1))/2; 
                    hold on                    
                    
                    
                    switch smoothHist
                        
                        case 'spline'

                            n = length(edges);
                            w = edges(2)-edges(1);
                            t = linspace(edges(1),edges(end),n+1);
                            dt = diff(t);
                            Fvals = cumsum([0,vals.*dt]);                
                            F = spline(t, [0, Fvals, 0]);
                            DF = fnder(F);
                            pts = fnplt(DF, 'r', 2);
                            set(gca, 'colororderindex', comp)
                            area(pts(1, :), pts(2, :), 'LineStyle', 'none', 'FaceAlpha', 1 / numComp)    
                            set(gca, 'colororderindex', comp)
                            plot(pts(1, :), pts(2, :), 'LineWidth', linewidth)  
                            
                        case 'gaussian'
                        
                            % fit gaussian 
                            try
                                [f, ~] = fit(edges', vals', 'gauss1');
                                smoothValid = true;
                            catch ERR
                                smoothValid = false;
                            end
                            % plot bars of actual distribution
%                             [vals, edges, ~] = histcounts(m, round(size(m, 1) / 4),...
%                                 'Normalization', 'probability');        
%                             set(gca, 'colororderindex', comp)
%                             bar(edges, vals, 1, 'EdgeColor', 'none',...
%                                 'FaceAlpha', .4);    
%                             hold on
                            % plot smoothed line
                            if smoothValid
                                lotsOfEdges = linspace(edges(1), edges(end), 1000);
                                set(gca, 'colororderindex', comp)
                                area(lotsOfEdges, f(lotsOfEdges), 'LineStyle', 'none', 'FaceAlpha', 1 / numComp)   
                                set(gca, 'colororderindex', comp)
                                plot(lotsOfEdges, f(lotsOfEdges), 'LineWidth', linewidth)
                            end
                            
                        case 'none'

                            % plot
                            set(gca, 'colororderindex', comp)
                            bar(edges, vals, 1, 'EdgeColor', 'none',...
                                'FaceAlpha', .4);                               
                            set(gca, 'colororderindex', comp)
                            plot(edges, vals, 'LineWidth', linewidth)
                            
                    end
                    
                end
                
                if showLegend && row == 1 && col == 1 
                    if ~smoothHist
                        legend(comp_u)
                    else
                        leg = cell(numComp * 2, 1);
                        compCounter = 1;
                        for legComp = 1:numComp
                            leg{compCounter} = sprintf('%s (raw)', comp_u{legComp});
                            compCounter = compCounter + 1;
                            leg{compCounter} = sprintf('%s (smoothed)', comp_u{legComp});
                            compCounter = compCounter + 1;
                        end
                        legend(leg, 'Location', 'best')
                    end
                end
                
                % settings
                xlabel(yLabel)
                ylabel('Probability')
                set(gca, 'YTick', [])
                ax = gca;
                ax.XRuler.Axle.LineWidth = axLineWidth;
                ax.YRuler.Axle.LineWidth = axLineWidth;  
                
                spIdx = spIdx + 1;
                
            end
        end
        
        % update ylim 
        if equalYLim
            for s = 1:2:spIdx - 1
                ylim(sp(s), [yl_min, yl_max])
            end
        end
        
        set(fig_barhist, 'Color', [1, 1, 1])
        set(fig_barhist, 'visible', 'on')     
        
    end        

    %% stacked bar/area chart

    if plotStacked && numComp >= 2
        
        fig_stack = figure('name', name);
        set(fig_stack, 'defaultaxesfontsize', fontsize);
        
        % loop through each comparison/row/col 
        spIdx = 1;
        for row = 1:numRow
            for col = 1:numCol
                
                subplot(numRow, numCol, spIdx);
                
%                 for comp = 1:numComp
                
                    % filter table for current row and col
                    idx_filt = row_s == row & col_s == col;% & comp_s == comp;
                    tmp = tab(idx_filt, :);

                    % make subs for comparison
                    [aoi_u, ~, aoi_s] = unique(tmp.aoi);
                    [aoiComp_u, ~, aoiComp_s] = unique(tmp.(compare));
                    numAOI = length(aoi_u);

                    % calculate mean and SEM for each comparison
                    bar_mu = accumarray([aoiComp_s, aoi_s],...
                        tmp.(measure), [], @nanmean);

                    % make title
                    if strcmpi(col_u{col}, 'NONE')
                        titleCol = '';
                    else
                        titleCol = sprintf('%s: %s', cols, upper(col_u{col}));
                    end
                    if strcmpi(row_u{row}, 'NONE')
                        titleRow = '';
                    else
                        titleRow = sprintf('%s: %s', rows, upper(row_u{row}));
                    end
                    titleStr = sprintf('%s %s', titleRow, titleCol);

                    % plot
                    h_sbar = bar(bar_mu, 0.5, 'stacked');
                    hold on
                    set(gca, 'ColorOrderIndex', 1)
                    x = repmat([1.25; 1.75], 1, numAOI);
                    area(x, bar_mu, 'FaceAlpha', 0.5, 'LineStyle', 'none')
                    xlim([0, numComp + 1])
    %                 h_bar.FaceColor = 'flat';
    %                 colours = lines(histNumComp);
    %                 for hc = 1:histNumComp
    %                     h_bar.CData(hc, :) = colours(hc, :);
    %                 end
    %                 h_bar.LineWidth = linewidth;
    %                 set(h_errbar, 'LineWidth', linewidth);

                    % settings
                    xlabel(xLabel)
                    ylabel(yLabel)
                    set(gca, 'XTickLabel', comp_u)
                    title(titleStr)
                    ax = gca;
                    ax.XRuler.Axle.LineWidth = axLineWidth;
                    ax.YRuler.Axle.LineWidth = axLineWidth; 
                    yl = [min(bar_mu - bar_sem), max(bar_mu + bar_sem)];
                    ylim([yl(1) - abs(yl(1) * .2), yl(2) + abs(yl(2) * .2)])
                    
                spIdx = spIdx + 1;
                
            end
        end
        set(fig_stack, 'Color', [1, 1, 1])        
        
        
        
        
    elseif plotStacked && numComp < 2
        
        warning('Cannot plot stacked bar chart unless comparing at least two variables.')

    end
       
    %% final settings

    set(gcf, 'Color', [1, 1, 1])
    hold off
    tilefigs
    
end

function [u, s, num] = prepareVariable(tab, var)

    if isempty(var)
        
        % if var is empty, return 'NONE' as a single unique value, and
        % vector of ones as the subscripts
        u = {'NONE'};
        numData = size(tab, 1);
        s = ones(numData, 1);
        num = 1;
        
    else
        
        % find unique vals of var, and subscripts, and do not sort (so that
        % pre-sorted tables can retain their ordering)
        [u, ~, s] = unique(tab.(var), 'stable');
        num = length(u);
        
        % if numeric, convert to cellstr
        if isnumeric(u) || islogical(u)
            u = arrayfun(@num2str, u, 'UniformOutput', false);
        end
        
        % convert any individual values that are numeric to string
        idx_isnum = cellfun(@isnumeric, u);
        u(idx_isnum) = cellfun(@num2str, u(idx_isnum),...
            'UniformOutput', false);  
    end
    
    % don't think we need this?
%     if isnumeric(u), u = num2cell(u); end
    
end

%% unused

%     %% t-test
%     
%     if doTTest && numComp ~= 2    
%         warning('t-test not performed since num comparions ~= 2.')
%     end
%         
%     res_t.lat_t = nan(numRow, numCol);
%     res_t.lat_p = nan(numRow, numCol);
%     res_t.lat_ci = nan(2, numRow, numCol);
%     res_t.lat_df = nan(numRow, numCol);
%     res_t.amp_t = nan(numRow, numCol);
%     res_t.amp_p = nan(numRow, numCol);
%     res_t.amp_ci = nan(2, numRow, numCol);
%     res_t.amp_df = nan(numRow, numCol);
% 
%     for row = 1:numRow
%         for col = 1:numCol
% 
%             % amp
%             [~, p, ci, stats] = ttest2(...
%                 meanamp{1, row, col},...
%                 meanamp{2, row, col});
%             res_t.amp_t(row, col) = stats.tstat;
%             res_t.amp_p(row, col) = p;
%             res_t.amp_ci(:, row, col) = ci;
%             res_t.amp_df(row, col) = stats.df;
% 
%             % lat
%             [~, p, ci, stats] = ttest2(...
%                 lat{1, row, col},...
%                 lat{2, row, col});
%             res_t.lat_t(row, col) = stats.tstat;
%             res_t.lat_p(row, col) = p;
%             res_t.lat_ci(:, row, col) = ci;
%             res_t.lat_df(row, col) = stats.df;
% 
%         end 
%     end
%             
%     %% make table
%     
%     if doTable && numComp == 2
%       
%         % make col labels
%         numElements = numRow * numCol;
%         colLab = cell(numElements, 1);
%         i = 1;
%         
%         tab_amp_mudiff = nan(numElements, 1);
%         tab_amp_sd = nan(numElements, 1);
%         tab_amp_se = nan(numElements, 1);
%         tab_amp_df1 = nan(numElements, 1);
%         tab_amp_t = nan(numElements, 1);
%         tab_amp_p = nan(numElements, 1);
%         tab_amp_d = nan(numElements, 1);
%         
%         tab_lat_mudiff = nan(numElements, 1);
%         tab_lat_sd = nan(numElements, 1);
%         tab_lat_se = nan(numElements, 1);
%         tab_lat_df1 = nan(numElements, 1);
%         tab_lat_t = nan(numElements, 1);
%         tab_lat_p = nan(numElements, 1);
%         tab_lat_d = nan(numElements, 1);
%         
%         for r = 1:numRow
%             for c = 1:numCol
%                 
%                 % column label
%                 colLab{i} = sprintf('%s_%s', row_u{r}, col_u{c});
%                 
%                 % amp
%                 N                   = sum(cellfun(@length, meanamp(:, r, c)));
%                 tab_amp_mudiff(i)   = nanmean(meanamp{1, r, c}) - nanmean(meanamp{2, r, c});
%                 tab_amp_sd(i)       = nanstd([cell2mat(meanamp(1, r, c)); cell2mat(meanamp(2, r, c))]);
%                 tab_amp_se(i)       = tab_amp_sd(i) / sqrt(N);
%                 tab_amp_df1(i)      = res_t.amp_df(r, c);
%                 tab_amp_t(i)        = res_t.amp_t(r, c);  
%                 tab_amp_p(i)        = res_t.amp_p(r, c);
%                 tab_amp_d(i)        = abs(tab_amp_mudiff(i) / tab_amp_sd(i));                 
%                 
%                 % lat
%                 N                   = sum(cellfun(@length, lat(:, r, c)));
%                 tab_lat_mudiff(i)   = nanmean(lat{1, r, c}) - nanmean(lat{2, r, c});
%                 tab_lat_sd(i)       = nanstd([cell2mat(lat(1, r, c)); cell2mat(lat(2, r, c))]);  
%                 tab_lat_se(i)       = tab_lat_sd(i) / sqrt(N);
%                 tab_lat_df1(i)      = res_t.lat_df(r, c);
%                 tab_lat_t(i)        = res_t.lat_t(r, c);  
%                 tab_lat_p(i)        = res_t.lat_p(r, c);
%                 tab_lat_d(i)        = abs(tab_lat_mudiff(i) / tab_lat_sd(i)); 
%                 
%                 i = i + 1;
%                         
%             end
%         end
%         
%         % make row labels
%         rowLab = {'comp', 'mean_diff', 'sd_pooled', 'se_pooled',...
%             'df', 't', 'p', 'd'};
%         
%         % make matrices of numeric data
%         tabData_amp = [tab_amp_mudiff, tab_amp_sd, tab_amp_se,...
%             tab_amp_df1, tab_amp_t, tab_amp_p, tab_amp_d];
%         tabData_lat = [tab_lat_mudiff, tab_lat_sd, tab_lat_se,...
%             tab_lat_df1, tab_lat_t, tab_lat_p, tab_lat_d];
%         
%         % make amp table
%         tab_stats_amp = array2table(tabData_amp, 'VariableNames',...
%             rowLab(2:end));
%         
%         % add column for comparison
%         tab_stats_amp.comp = colLab;
%         
%         % add column for measure
%         tab_stats_amp.measure = repmat({'amplitude'}, numElements, 1);
%         
%         % add column for component
%         tab_stats_amp.component = repmat({component}, numElements, 1);
%         
%         % add column for analysis name
%         tab_stats_amp.name = repmat({name}, numElements, 1);
%         
%         % reorder columns
%         tab_stats_amp = tab_stats_amp(:, [11, 10, 9, 8, 1:7]);
%         
%         % do the same for latency
%         tab_stats_lat = array2table(tabData_lat, 'VariableNames',...
%             rowLab(2:end));
%         tab_stats_lat.comp = colLab;
%         tab_stats_lat.measure = repmat({'latency'}, numElements, 1);
%         tab_stats_lat.component = repmat({component}, numElements, 1);
%         tab_stats_lat.name = repmat({name}, numElements, 1);
%         tab_stats_lat = tab_stats_lat(:, [11, 10, 9, 8, 1:7]);
%         
%         % cat tables
%         tab_stats = [tab_stats_amp; tab_stats_lat];
%         
%     elseif doTable && numComp ~= 2
%         
%         warning('Cannot produce summary table when number of comparisons > 2.')
%         tab_stats = table;
%         
%     else
%         
%         tab_stats = table;
%         
%     end    
%     
%     %% erps
%     
%     % make time vector for area plots
%     time_area = [time, fliplr(time)];
% 
%     if isempty(fig)
%         fig = figure('name', name, 'defaultaxesfontsize', fontsize);
%     else
%         clf
%         set(fig, 'defaultaxesfontsize', fontsize);
%     end
% 
%     % loop through each comparison/row/col 
%     spIdx = 1;
%     titleStr = cell(numRow, numCol);
%     for row = 1:numRow
%         for col = 1:numCol
%             
%             subplot(numRow, numCol, spIdx);
%             spIdx = spIdx + 1;
%             
%             arCol = feval(colMap, numComp);
%             
%             % legend
%             if isnumeric(comp_u)
%                 uCompStr = num2str(comp_u);
%             else
%                 uCompStr = comp_u;
%             end
%             for arComp = 1:numComp
%                 hold on
%                 dta = avg(:, arComp, row, col);
%                 
%                 if doDeTrend
%                     dta = detrend(dta);
%                 end
%                 
%                 plot(time, dta, 'linewidth', linewidth,...
%                         'color', arCol(arComp, :))
%             end         
%             if numComp > 1 && plotSEM && showLegend, legend(uCompStr, 'Location', 'NorthEast'), end
%             if numComp > 1 && ~plotSEM && showLegend, legend(uCompStr, 'Location', 'NorthEast'), end
%                 
%             for arComp = 1:numComp
%                 hold on
%                 if plotSEM
%                     
%                     dta = sem(:, arComp, row, col);
%                     
%                     if doDeTrend
%                         dta = detrend(dta);
%                     end
%                 
%                     ar = fill(time_area, dta, arCol(arComp, :));
%                     ar.FaceAlpha = SEMAlpha;
%                     ar.LineStyle = 'none';
%                 end
%             end
% 
%             plERP = gca;
%             
%             % title labels
%             str = '';
%             if numRow > 1, str = [str, ' | ', row_u{row}]; end
%             if numCol > 1, str = [str, ' | ', col_u{col}]; end
%             if strcmpi(str, '') || strcmpi(str, ' | ')
%                 str = name;
%             else 
%                 str = str(4:end);
%             end
%             
%             title(str)
%             titleStr{row, col} = str;
%             
%             % axis details
%             text(0.5, 0.025, 'Shaded areas +/- 2 SEM', 'units',...
%                 'normalized', 'horizontalalignment', 'center',...
%                 'verticalalignment', 'bottom', 'color', [.4, .4, .4])
%             xlabel('Time (s)')
%             ylabel('Amplitude (uV)')
%             set(gca, 'xgrid', 'on')
%             set(gca, 'xminorgrid', 'on')
%             set(gca, 'ygrid', 'on')
%             set(gca, 'ylim', [min(avg(:)) - 1, max(avg(:)) + 1]);
%             ax = gca;
%             ax.XRuler.Axle.LineWidth = axLineWidth;
%             ax.YRuler.Axle.LineWidth = axLineWidth;
%             
%         end
%     end
%     
%     set(gcf, 'Color', [1, 1, 1])
%     
%     %% boxplots
%     if plotBoxPlot
%         
%         % mean amp
%         fig_bp_amp = figure('name', [name, '_boxplot_meanamp']);
%         set(fig_bp_amp, 'defaultaxesfontsize', fontsize);
%         % ylim
%         ymin = min(vertcat(meanamp{:}));
%         ymax = max(vertcat(meanamp{:}));
%         yrange = ymax - ymin;
%         yl = [ymin, ymax];
%         % loop through each comparison/row/col 
%         spIdx = 1;
%         for row = 1:numRow
%             for col = 1:numCol
%                 % plot
%                 subplot(numRow, numCol, spIdx);
%                 hold on
%                 spIdx = spIdx + 1;
%                 for comp = 1:numComp
% %                     meanamp_bp = cell2mat(meanamp(:, row, col));
% %                     boxplot(meanamp_bp, comp_s, 'parent', fig_bp_amp);
%                     notBoxPlot(cell2mat(meanamp(comp, row, col)), comp, 'jitter', .5)
%                 end
%                 
%                 % append ttest results to title (if requested)
%                 if doTTest && ~isempty(res_t)
%                     str = sprintf('%s\nt(%d)=%.3f, p=%.3f',...
%                         titleStr{row, col},...
%                         res_t.amp_df(row, col),...
%                         res_t.amp_t(row, col),...
%                         res_t.amp_p(row, col));
%                 else
%                     str = titleStr{row, col};
%                 end
%                 
%                 % settings
%                 set(gca, 'xtick', 1:numComp)
%                 set(gca, 'xticklabel', comp_u)
%                 xlabel(compare, 'Interpreter', 'none')
%                 ylabel('Mean Amplitude (uV)')
%                 title(str)
%                 set(gca, 'ylim', yl);
%                 ax = gca;
%                 ax.XRuler.Axle.LineWidth = axLineWidth;
%                 ax.YRuler.Axle.LineWidth = axLineWidth;                
%                 
%             end
%         end
%         set(fig_bp_amp, 'Color', [1, 1, 1])
%         
%         % latency
%         fig_bp_lat = figure('name', [name, '_boxplot_latency']);
%         set(fig_bp_lat, 'defaultaxesfontsize', fontsize);
%         % loop through each comparison/row/col 
%         spIdx = 1;
%         for row = 1:numRow
%             for col = 1:numCol
%                 % plot
%                 sp = subplot(numRow, numCol, spIdx);
%                 hold on
%                 spIdx = spIdx + 1;
%                 for comp = 1:numComp
% %                     lat_bp = cell2mat(lat(:, row, col));
% %                     boxplot(lat_bp, comp_s, 'parent', sp);
%                     notBoxPlot(cell2mat(lat(comp, row, col)), comp, 'jitter', .5)
%                 end
%                 
%                 % append ttest results to title (if requested)
%                 if doTTest && ~isempty(res_t)
%                     str = sprintf('%s\nt(%d)=%.3f, p=%.3f',...
%                         titleStr{row, col},...
%                         res_t.lat_df(row, col),...
%                         res_t.lat_t(row, col),...
%                         res_t.lat_p(row, col));
%                 else
%                     str = titleStr{row, col};
%                 end
%                 
%                 % settings
%                 set(gca, 'xtick', 1:numComp)
%                 set(gca, 'xticklabel', comp_u)
%                 xlabel(compare, 'Interpreter', 'none')
%                 ylabel('Latency (s)')
%                 title(str)
%                 set(gca, 'ylim', [min(vertcat(lat{:})), max(vertcat(lat{:}))]);
%                 ax = gca;
%                 ax.XRuler.Axle.LineWidth = axLineWidth;
%                 ax.YRuler.Axle.LineWidth = axLineWidth;                
%                 
%             end
%         end
%         set(fig_bp_lat, 'Color', [1, 1, 1])
%         
%     end


%     %% back to back histograms
%     if plotB2BHist
%         
%      % mean amp
%         fig_b2bhist = figure('name', [name, '_boxplot_meanamp']);
%         set(fig_b2bhist, 'defaultaxesfontsize', fontsize);
%         % loop through each comparison/row/col 
%         spIdx = 1;
%         for row = 1:numRow
%             for col = 1:numCol
%                 % plot
%                 subplot(numRow, numCol, spIdx);
%                 hold on
%                 spIdx = spIdx + 1;
%                 b2bhist(cell2mat(meanamp(1, row, col)'), cell2mat(meanamp(2, row, col)'));
%                 % settings
%                 xlabel('Mean Amplitude (uV)')
%             end
%         end
%         set(fig_b2bhist, 'Color', [1, 1, 1])
%         
%         % latency
%         fig_hist2 = figure('name', [name, '_boxplot_latency']);
%         set(fig_hist2, 'defaultaxesfontsize', fontsize);
%         % loop through each comparison/row/col 
%         spIdx = 1;
%         for row = 1:numRow
%             for col = 1:numCol
%                 % plot
%                 subplot(numRow, numCol, spIdx);
%                 hold on
%                 spIdx = spIdx + 1;
%                 b2bhist(cell2mat(lat(1, row, col)'), cell2mat(lat(2, row, col)'));
%                 % settings
%                 xlabel('Latency (s)')
%             end
%         end
%         set(fig_hist2, 'Color', [1, 1, 1])        
%         
%     end
