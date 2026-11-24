# MATLAB 预处理对照

参考函数均在本目录维护一份，函数注释注明固定上游提交 `2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c`、文件位置与封装范围。

| 函数 | 出处与用途 |
| --- | --- |
| `original_preprocess(BOLD, TR)` | 上游 `StateSpace_modelling.m` 第 70–76、87–88 行：原样保留滤波设置和预处理语句，只添加函数参数、脑区循环及输出。 |
| `demean(x, dim)` | 上游 `functions/demean.m`，仅增加来源注释，用于去均值。 |
| `testPreprocessing` | 我们的测试驱动：构造数据，直接调用上述 MATLAB 函数，与安装包 `weida_preprocessing.preprocess_scan` 比较。 |

MATLAB R2025b + Signal Processing Toolbox 独立计算：去线性趋势 → 去均值 → 二阶 Butterworth 0.04–0.07 Hz 带通 → `filtfilt`。78 脑区、13/48/210/312 帧、TR 为 2.1/1.5/2.1/2.0 秒；包含所有首尾时间点，检查形状、有限值及输入不被修改，逐元素绝对容差 `1e-8`。

不执行 Hilbert、特征分解或聚类，不使用患者数据。不另存整份原脚本。

在仓库根目录执行：

```sh
uv sync --locked --python 3.11
matlab -batch "results = runtests('tests/matlab'); disp(table(results)); assertSuccess(results)"
```

## PhiID 读写

- `original_read_bold`：作者第 25–28 行原语句，仅加参数与输出；对照 78 脑区、1/13/210 帧及带/不带 BOM，绝对容差 1e-12。
- `write_table`：工程扩展，MATLAB `readtable` 验证字段、Unicode/引号、双精度数值、NaN/Inf。
- ID 配对、坏文件检查由 Python 单测覆盖；作者按文件顺序固定人数，无对应函数，不宣称 MATLAB 等价。
