function t = etTimeBuffer2Secs(tb)
%ETTIMEBUFFER2SECS Convert Task Engine microsecond timestamps to zeroed seconds.
if isempty(tb)
    t = [];
    return
end
t = double(tb(:, 1) - tb(1, 1)) / 1e6;
end
