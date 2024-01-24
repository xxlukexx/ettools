function avg_ind = etAverageAOIEvents(events, prePost, in, gaze, def)

    if ~exist('prePost', 'var') || isempty(prePost)
        prePost = [-1.0, 5.0];
    end
    
    numSubs = gaze.NumSubjects;
    numAOIs = size(in, 3);
    
    % segment all events

        numEvents = size(events, 1);
        seg = cell(numEvents, 1);
        t_seg = cell(numEvents, 1);
        for e = 1:numEvents

            % find timestamps of segment edges
            t1 = events{e, 2} + prePost(1);
            t2 = events{e, 2} + prePost(2);

            % convert to samples
            s1 = find(gaze.Time >= t1, 1);
            s2 = find(gaze.Time >= t2, 1);

            % segment AOI data, and time
            seg{e} = in(s1:s2, :, :);
            t_seg{e} = gaze.Time(s1:s2) - events{e, 2};

        end  
    
    % tabulate
    
        % find length of maximum segment
        lens = cellfun(@length, t_seg);
        mxLen = max(lens);
        
        % preallocate blank matrix for segments
        m = nan(mxLen, numSubs, numAOIs, numEvents);
        
        % use timestamps from longest segment as template
        t_m = t_seg{find(lens == mxLen, 1)};
    
        for e = 1:numEvents
            
            % refer to master timestamps to find sample edges for this
            % segment
            delta1 = abs(t_m - t_seg{e}(1));
            s1 = find(delta1 == min(delta1), 1);
            delta2 = abs(t_m - t_seg{e}(end));
            s2 = find(delta2 == min(delta2), 1);
            
            if s2 - s1 + 1 < size(seg{e}, 1)
                s1 = s1 - 1;
            end
            
            % place segment
            m(s1:s2, :, :, e) = seg{e};            
            
        end

    
    % group events by type
    
        [ev_u, ~, ev_s] = unique(events(:, 1));
        numEvTypes = length(ev_u);

        avg_ind = cell(numEvTypes, 1);
        avg_grand = cell(numEvTypes, 1);
        for et = 1:numEvTypes

            % filter by event type
            idx = ev_s == et;
            
            % ind avg
            avg_ind{et} = nanmean(m(:, :, :, idx), 4);
            
            % baseline ind avg
            bl_mu = nanmean(avg_ind{et}(t_m < 0, :, :, et));
            avg_ind{et} = avg_ind{et} - bl_mu;
            
            avg_grand{et} = shiftdim(nanmean(avg_ind{et}, 2), 2)';

        end

        figure
        plot(t_m, avg_grand{1})
        legend(def(:, 1))
    
end

   
%     % get timestamps of events
%     t_ev = cell2mat(events(:, 2));
%     
%     % adjust timestamps for pre/post
%     t1_ev = t_ev + prePost(1);
%     t2_ev = t_ev + prePost(2);
%     
%     % convert timestamps to samples
%     s_ev = arrayfun(@(x) find(gaze.Time >= x, 1), t_ev);
%     s1_ev = arrayfun(@(x) find(gaze.Time >= x, 1), t1_ev);
%     s2_ev = arrayfun(@(x) find(gaze.Time >= x, 1), t2_ev);
% %     numSamps = s2_ev - s1_ev + 1;
%     
%     % segment
%     seg = arrayfun(@(s1, s2) in(s1:s2, :, :), s1_ev, s2_ev,...
%         'UniformOutput', false);
%     t_seg = arrayfun(@(s1, s2, t_ev) gaze.Time(s1:s2) - t_ev, s1_ev, s2_ev, t_ev,...
%         'UniformOutput', false);