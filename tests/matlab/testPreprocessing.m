function tests = testPreprocessing
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
upstream = fullfile(root, 'tests', 'matlab', 'upstream', 'WEIDA-DBS');
source = fileread(fullfile(upstream, 'StateSpace_modelling.m'));
% Execute the author's statements verbatim, stopping before Hilbert/WEiDA.
settings = regexp(source, '(?m)^fnq=1/\(2\*TR\);[\s\S]*?^clear fnq flp fhi Wn k', 'match');
amplitude = regexp(source, '(?m)^\s*ts=demean\(detrend\(BOLD\(seed,:\)\)\);\s*signal_filt =filtfilt\(bfilt,afilt,ts\);', 'match');
assert(numel(settings)==1 && numel(amplitude)==1, ...
    'Expected exactly one original filter setup and amplitude block');
tc.TestData.settings = settings{1};
tc.TestData.amplitude = amplitude{1};
tc.TestData.oldPath = path;
addpath(fullfile(upstream, 'functions'), '-begin');
assert(strcmp(which('demean'), fullfile(upstream, 'functions', 'demean.m')));
if ispc
    interpreter = fullfile(root, '.venv', 'Scripts', 'python.exe');
else
    interpreter = fullfile(root, '.venv', 'bin', 'python');
end
pyenv('Version', interpreter, 'ExecutionMode', 'OutOfProcess');
insert(py.sys.path, int32(0), root);
% Test the installed package consumed by WEIDA and PhiID.
tc.TestData.processing = py.importlib.import_module('weida_preprocessing');
end

function teardownOnce(tc)
path(tc.TestData.oldPath);
end

function testOriginalAmplitudePipeline(tc)
rng(284, 'twister');
frames = [13 48 210 312];
trs = [2.1 1.5 2.1 2.0];
for c=1:numel(frames)
    TR = trs(c);
    t = (0:frames(c)-1)*TR;
    BOLD = 3+0.02*t+sin(2*pi*.055*t+(1:78)')+0.2*randn(78,frames(c));
    eval(tc.TestData.settings);
    expected = zeros(size(BOLD));
    for seed=1:size(BOLD,1)
        eval(tc.TestData.amplitude);
        expected(seed,:) = signal_filt;
    end
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
