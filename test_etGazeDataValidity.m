function tests = test_etGazeDataValidity
%TEST_ETGAZEDATAVALIDITY Analysis gaze excludes invalid/off-screen eyes.
tests = functiontests(localfunctions);
end

function testBinocularImportMasksEachEyeBeforeAveraging(testCase)
lx = [.2; -.1; .4; .4; .4; -.1];
ly = [.3;  .3; 1.2; .4; .4;  .3];
rx = [.6;  .6; .6; -2; .6; 1.1];
ry = [.7;  .7; .7; .7; .7;  .7];
leftMissing = [false; false; false; false; true; false];
rightMissing = false(6, 1);

gaze = etGazeDataBino;
gaze.Import(lx, ly, rx, ry, (0:5)', leftMissing, rightMissing);

verifyEqual(testCase, gaze.LeftMissing, ...
    [false; true; true; false; true; true])
verifyEqual(testCase, gaze.RightMissing, ...
    [false; false; false; true; false; true])
verifyEqual(testCase, gaze.Missing, ...
    [false; false; false; false; false; true])
verifyEqual(testCase, gaze.LeftX, [.2; nan; nan; .4; nan; nan])
verifyEqual(testCase, gaze.LeftY, [.3; nan; nan; .4; nan; nan])
verifyEqual(testCase, gaze.RightX, [.6; .6; .6; nan; .6; nan])
verifyEqual(testCase, gaze.RightY, [.7; .7; .7; nan; .7; nan])
verifyEqual(testCase, gaze.X, [.4; .6; .6; .4; .6; nan], ...
    'AbsTol', 1e-12)
verifyEqual(testCase, gaze.Y, [.5; .7; .7; .4; .7; nan], ...
    'AbsTol', 1e-12)
end

function testPreserveOptionKeepsPerEyeEvidenceButNotInvalidAverage(testCase)
gaze = etGazeDataBino('PreserveInvalidCoordinates', true);
gaze.Import(-1, .2, .5, .6, 0, false, false);

verifyEqual(testCase, gaze.LeftX, -1)
verifyEqual(testCase, gaze.LeftY, .2)
verifyTrue(testCase, gaze.LeftMissing)
verifyFalse(testCase, gaze.RightMissing)
verifyEqual(testCase, gaze.X, .5)
verifyEqual(testCase, gaze.Y, .6)
end

function testMonocularImportTreatsOnlyOutsideClosedUnitIntervalAsMissing(testCase)
x = [0; 1; -.01; 1.01; .5];
y = [1; 0; .5; .5; nan];
gaze = etGazeDataMono;
gaze.Import(x, y, (0:4)');

verifyEqual(testCase, gaze.Missing, [false; false; true; true; true])
verifyEqual(testCase, gaze.X, [0; 1; nan; nan; nan])
verifyEqual(testCase, gaze.Y, [1; 0; nan; nan; nan])
end

function testMonocularMainBufferAveragesOnlyValidEye(testCase)
mainBuffer = nan(2, 26);
mainBuffer(:, [13 26]) = 0;
mainBuffer(1, [7 8 20 21]) = [-1, -1, .6, .7];
mainBuffer(2, [7 8 20 21]) = [1, 1, .2, .3];
timeBuffer = uint64([0 0; 1e6 1e6]);

gaze = etGazeDataMono('mainBuffer', mainBuffer, 'timeBuffer', timeBuffer);

verifyEqual(testCase, gaze.X, [.6; .6], 'AbsTol', 1e-12)
verifyEqual(testCase, gaze.Y, [.7; .65], 'AbsTol', 1e-12)
verifyFalse(testCase, any(gaze.Missing))
end
