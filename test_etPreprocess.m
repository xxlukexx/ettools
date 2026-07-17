function tests = test_etPreprocess
tests = functiontests(localfunctions);
end

function testMissingAndOffscreenSamples(testCase)
buffer = zeros(4, 26);
buffer(:, [7, 8, 20, 21]) = .5;
buffer(1, [13, 26]) = 4;
buffer(2, 7) = 1.2;
buffer(3, 20) = -.1;

[actual, ~, info] = etPreprocess(buffer, ...
    'removeoffscreen', true, 'removemissing', true);

verifyTrue(testCase, all(isnan(actual(1, [1:12, 14:25]))))
verifyEqual(testCase, actual(1, [13, 26]), [4, 4])
verifyTrue(testCase, all(isnan(actual(2, 7:8))))
verifyEqual(testCase, actual(2, 13), 4)
verifyEqual(testCase, actual(2, 20:21), [.5, .5])
verifyTrue(testCase, all(isnan(actual(3, 20:21))))
verifyEqual(testCase, actual(3, 26), 4)
verifyEqual(testCase, info.propValL, .5)
verifyEqual(testCase, info.propValR, .5)
verifyEqual(testCase, info.propVal, .75)
end

function testNoOptionsPreservesBuffer(testCase)
buffer = rand(5, 26);
buffer(:, [13, 26]) = 0;
[actual, time, info] = etPreprocess(buffer);
verifyEqual(testCase, actual, buffer)
verifyEmpty(testCase, time)
verifyEqual(testCase, info.propVal, 1)
end

function testTimeBufferConversion(testCase)
timeBuffer = uint64([1000000, 0; 1005000, 0; 1010000, 0]);
verifyEqual(testCase, etTimeBuffer2Secs(timeBuffer), [0; .005; .01], ...
    'AbsTol', 1e-12)
end
