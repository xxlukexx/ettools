function tests = test_etResample3
%TEST_ETRESAMPLE3 Regular-bin gaze resampling behavior.

    tests = functiontests(localfunctions);
end

function testDownsamplingAveragesContinuousValues(testCase)

    mb = zeros(4, 26);
    mb(:, 1:12) = repmat([1; 3; 5; 7], 1, 12);
    mb(:, 14:25) = repmat([2; 4; 6; 8], 1, 12);
    mb(:, 13) = [0; 1; 0; 1];
    mb(:, 26) = [1; 0; 1; 0];
    tb = uint64([1000000, 1; 1010000, 2; 1020000, 3; 1030000, 4]);

    [actual, actualTimeBuffer, actualTime] = etResample3(mb, tb, 50);

    verifyEqual(testCase, actualTime, [0; .02], 'AbsTol', 1e-12)
    verifyEqual(testCase, actual(:, 1:12), repmat([2; 6], 1, 12))
    verifyEqual(testCase, actual(:, 14:25), repmat([3; 7], 1, 12))
    verifyEqual(testCase, actual(:, 13), [0; 0])
    verifyEqual(testCase, actual(:, 26), [0; 0])
    verifyEqual(testCase, actualTimeBuffer, tb([1, 3], :))
end

function testUpsamplingCopiesNearestSamplesIntoEmptyBins(testCase)

    mb = zeros(2, 26);
    mb(1, [1:12, 14:25]) = 1;
    mb(2, [1:12, 14:25]) = 2;
    tb = uint64([2000000, 11; 2100000, 12]);

    [actual, actualTimeBuffer, actualTime] = etResample3(mb, tb, 20);

    verifyEqual(testCase, actualTime, [0; .05; .1], 'AbsTol', 1e-12)
    verifyEqual(testCase, actual(:, 1), [1; 1; 2])
    verifyEqual(testCase, actual(:, 14), [1; 1; 2])
    verifyEqual(testCase, actualTimeBuffer, tb([1, 1, 2], :))
end

function testValidityAggregationIgnoresNaNs(testCase)

    mb = zeros(2, 26);
    mb(:, 13) = [nan; 2];
    mb(:, 26) = [3; nan];
    tb = uint64([3000000, 1; 3010000, 2]);

    actual = etResample3(mb, tb, 50);

    verifyEqual(testCase, actual(1, 13), 2)
    verifyEqual(testCase, actual(1, 26), 3)
end

function testExtraInputColumnsAreRemoved(testCase)

    mb = zeros(2, 30);
    tb = uint64([4000000, 1; 4010000, 2]);

    actual = etResample3(mb, tb, 50);

    verifySize(testCase, actual, [1, 26])
end

function testEmptyInputsReturnEmptyOutputs(testCase)

    [mb, tb, time] = etResample3([], [], 50);

    verifyEmpty(testCase, mb)
    verifyEmpty(testCase, tb)
    verifyEmpty(testCase, time)
end
