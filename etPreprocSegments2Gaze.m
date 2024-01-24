function gaze = etPreprocSegments2Gaze(data)

    for s = length(data.Segments):-1:1
        gaze(s) = etGazeDataBino('mainBuffer', data.Segments(s).MainBuffer,...
            'timeBuffer', data.Segments(s).TimeBuffer, 'eventBuffer',...
            data.Segments(s).EventBuffer);
    end

end