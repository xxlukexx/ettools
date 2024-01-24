function avg = etSummarisePupilSizeByLook(looks, baselineDurS, dataDurS)

    avg = nan(size(looks));
    for i = 1:numel(looks)
        
        % get indices of baseline, and averaging period
        baselineSamps = looks{i}.Time <= baselineDurS;
        bl = looks{i}.Pupil(baselineSamps);
        bl_mu = nanmean(bl);
        
        dataSamples = looks{i}.Time > baselineDurS &...
            looks{i}.Time <= (dataDurS + baselineDurS);
        data = looks{i}.Pupil(dataSamples);
        data_mu = nanmean(data);
        
        avg(i) = data_mu - bl_mu;
        
    end

end