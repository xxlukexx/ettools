numSubs = 100;
duration = 10;
fs = 120;
file_aoi = '/Users/luke/Google Drive/Experiments/popout/aoi_masks_te/POPOUT1_dilated.png';
load('/Users/luke/Google Drive/Experiments/popout/aoi_masks_te/aoi_def.mat')
numAOIs = size(def, 1);

% make sham data
fprintf('Making sham gaze data...\n')
[gaze, in, res_sham] = etMakeShamData(numSubs, duration, fs, file_aoi, def);

% score in visualiser    
fprintf('Scoring sham gaze data in visualiser...\n')
clear port stim et aoi
sca
port = teViewport;
port.LayerClickBox = port.Size;
port.LayerClickBox(4) = port.LayerClickBox(4) - 250;

tl = vpaTimeline_EyeTracking;
tl.Type = 'ui';
tl.ZPosition = 100;
port.Viewpane('timeline') = tl;

aoi = vpaStaticAOI;
aoi.ZPosition = 2;
aoi.Image = file_aoi;
port.Viewpane('aoi') = aoi;
aoi.AOIDefinition = def;
aoi.Alpha = .7;
tl.AOI = aoi;
aoi.InterpolateResultsSecs = 0;
aoi.TriggerToleranceSecs = 0;

et = vpaEyeTrackingData;
et.ZPosition = 3;
et.Import(gaze)
et.Timeline = tl;
tl.Duration = 10;
port.Viewpane('gaze') = et;
et.DrawHeatmap = false;
aoi.Score(et)

% calculate prop valid between etGazeData object and sham input
diff_val = gaze.PropValid - res_sham.propValid;

    fprintf('Test of etGazeData prop valid (sham created vs etGazeData calculated): ')
    if all(diff_val(:) == 0)
        fprintf('PASSED\n')
    else
        fprintf(2, 'FAILED\n')
    end

% calculate metrics for sham data (tests etSummariseAOI independently)
[smry_sham, tall_sham, wide_sham] = etSummariseAOI(in, gaze, def);
samps_sham = arrayfun(@(x) x.samplesInAOI, smry_sham)';
prop_sham = arrayfun(@(x) x.propInAOI, smry_sham)';

% calculate metrics for scored visualiser data (tests visualiser)
[smry_vis, tall_vis, wide_vis] = etSummariseAOI(aoi.In, gaze, def);
samps_vis = arrayfun(@(x) x.samplesInAOI, smry_vis)';
prop_vis = arrayfun(@(x) x.propInAOI, smry_vis)';

% compare calculated with sham
    diff_samps_sham = res_sham.sampsAOI - samps_sham;
    fprintf('Test of etSummariseAOI SAMPLES (sham created vs sham calculated): ')
    if all(diff_samps_sham(:) == 0)
        fprintf('PASSED\n')
    else
        fprintf(2, 'FAILED\n')
    end

    diff_prop_sham = res_sham.propAOI - prop_sham;
    fprintf('Test of etSummariseAOI PROP (sham created vs sham calculated): ')
    if all(diff_prop_sham(:) == 0)
        fprintf('PASSED\n')
    else
        fprintf(2, 'FAILED\n')
    end

% compare visualise with sham
    diff_samps_vis = res_sham.sampsAOI - samps_vis;
    fprintf('Test of visualiser SAMPLES (sham created vs visualiser calculated): ')
    if all(diff_samps_vis(:) == 0)
        fprintf('PASSED\n')
    else
        fprintf(2, 'FAILED\n')
    end
    
    diff_prop_vis = res_sham.propAOI - prop_vis;
    fprintf('Test of visualiser PROP (sham created vs visualiser calculated): ')
    if all(diff_prop_vis(:) == 0)
        fprintf('PASSED\n')
    else
        fprintf(2, 'FAILED\n')
    end    

clear port stim et aoi
sca
