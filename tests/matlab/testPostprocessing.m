function tests = testPostprocessing
% 固定算法输出后，与原作者后处理函数逐项比较，不运行 eigs 或聚类。
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
tc.TestData.model = py.importlib.import_module('weida');
tc.TestData.metrics = py.importlib.import_module('downstream');
end

function a = array(values)
a = py.numpy.array(py.json.loads(jsonencode(values)));
end

function a = integers(values)
a = py.numpy.atleast_1d(py.numpy.array(py.json.loads(jsonencode(values)), pyargs('dtype', 'int64')));
end

function values = numbers(a)
values = jsondecode(char(py.json.dumps(a.tolist())));
end

function verifyArray(tc, actual, expected)
% 单独比较缺失掩码，再比较有限值，避免 JSON 对 NaN 的非标准编码。
tc.verifyEqual(logical(numbers(py.numpy.isnan(actual))), isnan(expected));
expected(isnan(expected)) = 0;
tc.verifyEqual(numbers(py.numpy.nan_to_num(actual)), expected, 'AbsTol', 1e-12);
end

function testOccupancyAndDwellTime(tc)
sequences = {[1 1 2 2 3 3 2 2 1 1 4 4 1], [2 2 2 1 1 2], [1 1 1 1], [4]};
K = 4; TR = 2.1;
expectedP = zeros(numel(sequences),K); expectedLT = nan(numel(sequences),K);
for s=1:numel(sequences)
    [expectedP(s,:), expectedLT(s,:)] = original_state_occupancy(sequences{s}(:), K, TR);
end
labels = [sequences{:}]; lengths = cellfun(@numel,sequences);
r = tc.TestData.model.analyze_state_occupancy(integers(labels-1), integers(lengths), int32(K), TR);
verifyArray(tc, r{1}, expectedP);
verifyArray(tc, r{2}, expectedLT);
end

function testScanTransitionsAndFrequencyOrder(tc)
sequences = {[1 1 2 2 1 3 3], [2 2 2 2], [4]};
labels = [sequences{:}]; lengths = cellfun(@numel,sequences); K = 4;
P = zeros(3,K); LT = nan(3,K);
for s=1:3
    [P(s,:),LT(s,:)] = original_state_occupancy(sequences{s}(:),K,2.1);
end
r = tc.TestData.model.calculate_state_dynamics(integers(labels-1),array(eye(K)), ...
    integers(lengths),array(P),array(LT),py.list({'A','A','B'}));
% 原脚本第 4 节：按全体帧的状态占比降序排列。
ProbC = zeros(1,K);
for c=1:K
    ProbC(c)=mean(labels==c);
end
[~,ind_sort]=sort(ProbC,'descend');
tc.verifyEqual(numbers(r{'order'})+1,ind_sort(:));
for s=1:3
    IDX = sequences{s}(:);
    for c=1:K
        IDX(sequences{s}(:)==ind_sort(c))=c;
    end
    expected = original_transitions(IDX,K);
    verifyArray(tc, r{'transitions'}.take(int32(s-1),pyargs('axis',int32(0))),expected);
end
end

function testFixedCenterAssignment(tc)
rng(72,'twister');
centers = [0 0 0; 10 0 0; 0 10 0; 0 0 10];
features = [centers; centers+0.1*randn(4,3); 5 0 0; 0 5 0];
expected = original_predict_states(features,centers);
r = tc.TestData.metrics.predict_states(array(features),array(centers),integers([4 6]));
tc.verifyEqual(numbers(r{'states'})+1,expected);
tc.verifyEqual(expected(9:10),[1;1]); % 平距时选原始顺序中的第一个中心。
end

function testMetastability(tc)
rng(73,'twister');
features=randn(66,78); centers=randn(4,78); lengths=[37 29];
labels=mod((1:66)',4)+1; expected=zeros(2,4); offset=0;
for s=1:2
    idx=offset+(1:lengths(s)); offset=offset+lengths(s);
    expected(s,:)=original_metastability(features(idx,:),labels(idx),centers');
end
actual=tc.TestData.metrics.state_metastability(array(features),integers(labels-1),array(centers),integers(lengths));
verifyArray(tc,actual,expected);
end

function testUndefinedMetastabilityPolicy(tc)
% 原公式单帧状态除零、缺席状态返回零；Python 显式将两者视为不可估计。
features=[1 2; 3 4; 9 9]; labels=[1;1;2]; centers=zeros(3,2);
expected=original_metastability(features,labels,centers');
actual=tc.TestData.metrics.state_metastability(array(features),integers(labels-1),array(centers),integers(3));
tc.verifyEqual(expected(1),15);
tc.verifyTrue(isinf(expected(2)));
tc.verifyEqual(expected(3),0);
tc.verifyEqual(logical(numbers(py.numpy.isnan(actual))),[false true true]);
tc.verifyEqual(numbers(py.numpy.nan_to_num(actual)),[15 0 0]);
end
