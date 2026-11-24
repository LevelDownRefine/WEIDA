function tests = testPreprocessing
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ispc
    interpreter = fullfile(root, '.venv', 'Scripts', 'python.exe');
else
    interpreter = fullfile(root, '.venv', 'bin', 'python');
end
pyenv('Version', interpreter, 'ExecutionMode', 'OutOfProcess');
insert(py.sys.path, int32(0), root);
% 测试 WEIDA 和 PhiID 共用的已安装预处理包。
tc.TestData.processing = py.importlib.import_module('weida_preprocessing');
end

function testOriginalAmplitudePipeline(tc)
rng(284, 'twister');
frames = [13 48 210 312];
trs = [2.1 1.5 2.1 2.0];
for c=1:numel(frames)
    TR = trs(c);
    t = (0:frames(c)-1)*TR;
    BOLD = 3+0.02*t+sin(2*pi*.055*t+(1:78)')+0.2*randn(78,frames(c));
    expected = original_preprocess(BOLD, TR);
    input = py.numpy.array(py.json.loads(jsonencode(BOLD)));
    original = numbers(input);
    actual = numbers(tc.TestData.processing.preprocess_scan(input, pyargs( ...
        'remove_trend', true, 'demean', true, ...
        'band', py.tuple({.04, .07}), 'tr', TR)));
    tc.verifySize(actual, size(BOLD));
    tc.verifyTrue(all(isfinite(actual(:))));
    tc.verifyEqual(actual, expected, 'AbsTol', 1e-8);
    tc.verifyEqual(numbers(input), original);
    fprintf('Original WEIDA preprocessing: frames=%d TR=%g max_error=%.17g\n', ...
        frames(c), TR, max(abs(actual(:)-expected(:))));
end
end

function values = numbers(a)
values = jsondecode(char(py.json.dumps(a.tolist())));
end
