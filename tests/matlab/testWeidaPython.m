function tests = testWeidaPython
% MATLAB 直接调用项目的 Python 函数。每个测试的输入和预期值均在本文件中。
tests = functiontests(localfunctions);
end

function setupOnce(tc)
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ispc
    interpreter = fullfile(root, '.venv', 'Scripts', 'python.exe');
else
    interpreter = fullfile(root, '.venv', 'bin', 'python');
end
assert(isfile(interpreter), '请先在仓库根目录执行 uv sync --locked');
pyenv('Version', interpreter, 'ExecutionMode', 'OutOfProcess');
insert(py.sys.path, int32(0), root);
tc.TestData.model = py.importlib.import_module('weida');
tc.TestData.metrics = py.importlib.import_module('downstream');
end

function a = array(literal)
% 用内存中的嵌套列表显式保留矩阵形状与整数标签类型。
a = py.numpy.array(py.json.loads(literal));
end

function values = numbers(a)
values = jsondecode(char(py.json.dumps(a.tolist())));
end

function testRankOneFeatures(tc)
% 相同的四条非恒定信号对应全 1 相位同步矩阵，归一化后每个权重为 -1/2。
signal = sin((0:39)*0.7);
scan = array(jsonencode(repmat(signal,4,1)));
result = tc.TestData.model.compute_weighted_eigenvectors(py.list({scan}));
tc.verifyEqual(numbers(result{1}), -0.5*ones(40,4), 'AbsTol',1e-10);
end

function testKnownClusterCenters(tc)
x = array('[[0],[2],[10],[12]]');
r = tc.TestData.model.cluster_brain_states(x,pyargs('k',int32(2),'replicates',int32(8)));
tc.verifyEqual(sort(numbers(r{'centers'})), [1;11], 'AbsTol',1e-12);
tc.verifyEqual(numbers(r{'sumd'}), [2;2], 'AbsTol',1e-12);
end

function testCompleteRunsOnly(tc)
r = tc.TestData.model.analyze_state_occupancy(array('[0,0,1,1,0,1]'),array('[6]'),int32(2),2.0);
tc.verifyEqual(numbers(r{1}), [0.5,0.5], 'AbsTol',1e-12);
tc.verifyEqual(numbers(r{2}), [2,4], 'AbsTol',1e-12);
end

function testScanBoundaryDoesNotCreateTransition(tc)
labels = array('[0,0,1,1]'); lengths = array('[2,2]');
r = tc.TestData.model.calculate_state_dynamics(labels,array('[[1,0],[0,1]]'),lengths, ...
    array('[[1,0],[0,1]]'),array('[[0,0],[0,0]]'),py.list({'A','A'}));
groups = r{'group_transitions'};
tc.verifyEqual(numbers(groups{'A'}),eye(2));
end

function testFixedCentersKeepTheirOrder(tc)
r = tc.TestData.metrics.predict_states(array('[[10],[10],[0],[5]]'),array('[[0],[10]]'),array('[4]'));
tc.verifyEqual(numbers(r{'states'}),[1;1;0;0]);
end

function testMetastability(tc)
r = tc.TestData.metrics.state_metastability(array('[[1,2],[3,4]]'),array('[0,0]'),array('[[0,0]]'),array('[2]'));
tc.verifyEqual(double(r.item()),15); % (1+4+9+16)/(2*(2-1))
end

function testFunctionalSegregation(tc)
% 模块内相关均为 0.6，模块间相关全为零。使用四条正交的零均值序列。
h = hadamard(8); x = [h(2,:); .6*h(2,:)+.8*h(3,:); h(4,:); .6*h(4,:)+.8*h(5,:)];
r = tc.TestData.metrics.functional_segregation(array(jsonencode(x)),array('[1,1,-1,-1]'));
tc.verifyEqual(double(r{'segregation'}),atanh(.6),'AbsTol',1e-12);
tc.verifyEqual(numbers(r{'nodal'}),atanh(.6)*ones(4,1),'AbsTol',1e-12);
end

function testPairById(tc)
r = tc.TestData.metrics.paired_comparison(array('[4,0,1,0,2,0]'),py.list({'c','a','a','c','b','b'}),py.list({'ON','OFF','ON','OFF','ON','OFF'}));
tc.verifyEqual(numbers(r{'t'}),sqrt(7),'AbsTol',1e-12); % ON-OFF 分别为 1、2、4
end
