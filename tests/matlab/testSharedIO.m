function tests = testSharedIO
% 仅测试 PhiID 所需读写；不运行 WEiDA 状态算法。
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
tc.TestData.io = py.importlib.import_module('common.io');
end

function testOriginalCSVReading(tc)
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
for frames = [1 13 210]
    path = fullfile(folder, 'bold.csv');
    BOLD = reshape(sin(1:78*frames), 78, frames);
    for bom = [false true]
        fid = fopen(path, 'w', 'n', 'UTF-8');
        if bom, fwrite(fid, uint8([239 187 191]), 'uint8'); end
        fprintf(fid, 'index');
        fprintf(fid, ',ROI%d', 1:78); fprintf(fid, '\n');
        for t = 1:frames
            fprintf(fid, '%d', t+100);
            fprintf(fid, ',%.17g', BOLD(:,t)); fprintf(fid, '\n');
        end
        fclose(fid);
        expected = original_read_bold(path);
        result = tc.TestData.io.read_bold_csv(path);
        actual = double(result{1});
        tc.verifySize(actual, [78 frames]);
        tc.verifyEqual(actual, expected, 'AbsTol', 1e-12);
        tc.verifyEqual(actual, BOLD, 'AbsTol', 1e-12);
        tc.verifyEqual(int64(py.len(result{2})), int64(78));
    end
end
end

function testNamedTableReadByMATLAB(tc)
% 带表头、文字列是 PhiID 的工程格式；用 MATLAB readtable 检查互操作。
folder = tempname; mkdir(folder);
cleanup = onCleanup(@() rmdir(folder, 's'));
path = fullfile(folder, 'statistics.csv');
rows = py.list({ ...
    py.dict(pyargs('atom', '冗余,"信息"', 'value', pi, 'p', 1e-12)), ...
    py.dict(pyargs('atom', 'rtr', 'value', -2.5, 'p', py.float('nan'))), ...
    py.dict(pyargs('atom', 'rtx', 'value', py.float('inf'), 'p', 0.0))});
tc.TestData.io.write_table(path, rows);
t = readtable(path, 'Encoding', 'UTF-8', 'TextType', 'string');
tc.verifyEqual(height(t), 3);
tc.verifyEqual(t.Properties.VariableNames, {'atom', 'value', 'p'});
tc.verifyEqual(t.atom, ["冗余,""信息"""; "rtr"; "rtx"]);
tc.verifyEqual(t.value(1:2), [pi; -2.5], 'AbsTol', 1e-14);
tc.verifyEqual(t.value(3), Inf);
tc.verifyEqual(t.p([1 3]), [1e-12; 0], 'AbsTol', 1e-20);
tc.verifyTrue(isnan(t.p(2)));
end
