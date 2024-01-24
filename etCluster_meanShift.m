function out = etCluster_meanShift(gaze, bandwidth, vars, scale_vars)

    if ~exist('gaze', 'var') || isempty(gaze) || ~isa(gaze, 'etGazeData')
        error('Must supply ''gaze'', a teGazeData object.')
    end
    
    if ~exist('bandwidth', 'var') || isempty(bandwidth)
        bandwidth = 0.1;
        warning('''bandwidth'' was not specified, so has defaulted to 0.1.')
    end
    
    if ~exist('vars', 'var') 
        vars = [];
    end
    
    if ~exist('scale_vars', 'var')
        scale_vars = 1;
    end
    
    % reshape x, y
    m_x = reshape(gaze.X, 1, []);
    m_y = reshape(gaze.Y, 1, []);    
    
    % (optionally) format additional variables 
    if ~isempty(vars)
        if ~isvector(vars) || ~isnumeric(vars) || size(vars, 1) ~= gaze.NumSubjects
            error('Additional variables must be a vector or matrix with one row per subject and one each variable in a separate column.')
        end
        numVars = size(vars, 2);
        if numVars > 1
            error('More than one additional variable not supported.')
        end
        v = repmat(vars', gaze.NumSamples, 1);
        m_v = reshape(v, 1, []) .* scale_vars;
        m_gaze = [m_x; m_y; m_v];
    else
        numVars = 0;
        m_gaze = [m_x; m_y];
    end
       
    % remove missing
    missing                         = any(isnan(m_gaze), 1);
    m_gaze(:, missing)              = [];
    gNumPoints                      = size(m_gaze, 2);
    
    % cluster
    tic
    [cluster_centre, cluster_idx, cluster_members] =...
        HGMeanShiftCluster(m_gaze, bandwidth, 'gaussian');
    toc
    
    % validate clusters - must have > 1% of gaze samples
    clNumPoints                     = cellfun(@(x) size(x, 2), cluster_members);
    clPropPoints                    = clNumPoints ./ gNumPoints;
    clVal                           = clPropPoints >= 0.01;
    numClus                         = length(clVal);
    numValClus                      = sum(clVal);
    
    %store
    out.gaze                        = gaze;
    out.gaze_matrix                 = m_gaze;
    out.numGazePoints               = gNumPoints;
    out.numClusters                 = numClus;
    out.numValidClusters            = numValClus;
    out.cluster_centre              = cluster_centre;
    out.cluster_idx                 = cluster_idx;
    out.cluster_members             = cluster_members;
    out.cluster_validity            = clVal;
    out.cluster_propGaze            = clPropPoints;
    out.cluster_numGazePoints       = clNumPoints;
    
end