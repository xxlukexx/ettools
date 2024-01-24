function etGrandAverageAOI(avg_ind, groups)

    [grp_u, ~, grp_s] = unique(groups);
    numGroups = length(grp_u);
    
    numEventTypes = length(avg_ind);
    ga = cell(numGroups, numEventTypes);
    
    for g = 1:numGroups
        
        for e = 1:numEventTypes
    
            idx = grp_s == g;
            ga{g, e} = shiftdim(nanmean(avg_ind{e}(:, idx, :), 2), 2)';
            
        end
        
    end

    
    





end