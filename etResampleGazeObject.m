function gaze_r = etResampleGazeObject(gaze, fs_new)

    [mb, tb] = gaze.ExportTobiiAnalytics;
    [mb_r, tb_r] = etResample2(mb, tb, fs_new);
    gaze_r = etGazeDataBino('mainBuffer', mb_r, 'timeBuffer', tb_r);
    
end