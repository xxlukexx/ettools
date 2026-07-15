function tests = test_etGazeDataPupil

tests = functiontests(localfunctions);

end


function testSharedMissingMaskIsBackwardCompatible(testCase)

gaze = localGaze;
left = [3; 4; 5; 6];
right = [5; 6; 7; 8];
missing = [false; true; false; true];

gaze.ImportPupil(left, right, missing);

verifyEqual(testCase, gaze.LeftPupilMissing, missing)
verifyEqual(testCase, gaze.RightPupilMissing, missing)
verifyEqual(testCase, gaze.PupilMissing, missing)
verifyEqual(testCase, gaze.Pupil, [4; nan; 6; nan])

end


function testIndependentMissingMasksControlAverage(testCase)

gaze = localGaze;
left = [3; 4; 5; 6];
right = [5; 6; 7; 8];
leftMissing = [false; true; false; true];
rightMissing = [true; false; false; true];

gaze.ImportPupil(left, right, leftMissing, rightMissing);

verifyEqual(testCase, gaze.LeftPupilMissing, leftMissing)
verifyEqual(testCase, gaze.RightPupilMissing, rightMissing)
verifyEqual(testCase, gaze.PupilMissing, [false; false; false; true])
verifyEqual(testCase, gaze.Pupil, [3; 6; 6; nan])

end


function testTaskEngine2ExportPreservesPerEyeValidity(testCase)

gaze = localGaze;
leftMissing = [false; true; false; true];
rightMissing = [true; false; false; true];
gaze.ImportPupil([3; 4; 5; 6], [5; 6; 7; 8], ...
    leftMissing, rightMissing);

buffer = gaze.ExportTaskEngine2;

verifyEqual(testCase, logical(buffer(:, 6)), ~leftMissing)
verifyEqual(testCase, logical(buffer(:, 21)), ~rightMissing)
verifyEqual(testCase, buffer(:, 5), gaze.LeftPupil)
verifyEqual(testCase, buffer(:, 20), gaze.RightPupil)

end


function testSegmentPreservesPerEyeValidity(testCase)

gaze = localGaze;
leftMissing = [false; true; false; true];
rightMissing = [true; false; false; true];
gaze.ImportPupil([3; 4; 5; 6], [5; 6; 7; 8], ...
    leftMissing, rightMissing);

segment = gaze.SegmentBySample(2, 4);

verifyEqual(testCase, segment.LeftPupilMissing, leftMissing(2:4))
verifyEqual(testCase, segment.RightPupilMissing, rightMissing(2:4))
verifyEqual(testCase, segment.Pupil, [6; 6; nan])

end


function testLegacyStoredObjectFallsBackToSharedMask(testCase)

gaze = localGaze;
sharedMissing = [false; true; false; true];
gaze.LeftPupil = [3; 4; 5; 6];
gaze.RightPupil = [5; 6; 7; 8];
gaze.Pupil = [4; 5; 6; 7];
gaze.PupilMissing = sharedMissing;
gaze.LeftPupilMissing = [];
gaze.RightPupilMissing = [];

buffer = gaze.ExportTaskEngine2;
segment = gaze.SegmentBySample(2, 4);

verifyEqual(testCase, logical(buffer(:, 6)), ~sharedMissing)
verifyEqual(testCase, logical(buffer(:, 21)), ~sharedMissing)
verifyEqual(testCase, segment.LeftPupilMissing, sharedMissing(2:4))
verifyEqual(testCase, segment.RightPupilMissing, sharedMissing(2:4))

end


function gaze = localGaze

gaze = etGazeDataBino;
time = (0:3)' / 100;
gaze.Import(0.1 * ones(4, 1), 0.2 * ones(4, 1), ...
    0.3 * ones(4, 1), 0.4 * ones(4, 1), time, ...
    false(4, 1), false(4, 1), false(4, 1), time);

end
