# 单元测试、Octave 对照与 MATLAB 委托测试

所有命令在仓库根目录运行。测试使用合成数据，不需要上传真实受试者数据。

## 1. 单元测试（UT）

```bash
uv run --locked python -m unittest discover -s tests -t . -v
```

单元测试覆盖：完整状态段的持续时间、缺失状态和全程单一状态、扫描边界隔离、按扫描平均的组转移概率、状态排序一致性、非法长度与标签、全零及非有限输入、不同扫描长度、交接包生成和校验。

预期值来自可手算样例；测试不依赖 MATLAB 或 Octave。

## 2. 当前的 Octave 对照

```bash
uv run --locked python validation/compare.py --octave C:/path/to/octave-cli.exe --output validation_output
```

比较器生成固定种子的 78 脑区、210 时间点数据，调用 `reference.m`，输出中间 `.mat` 文件、日志和 `report.json`。

参考脚本按原公式独立实现，并非未经修改的整份上游脚本：数据读取改为固定测试输入；保留原脚本的求解器符号和最终向量定号顺序；聚类使用分离的模拟簇与共享初始中心，而非随机重启 300 次；状态统计使用预设标签，涵盖缺失状态、全程单一状态及首尾不完整段。

检查内容：

- 滤波系数、滤波信号和相位。
- 导出 MATLAB/Octave 原始特征对，调用生产函数 `weight_eigenvectors` 单独验证加权及最终定号，不在 Python 重新分解或匹配符号。
- 原符号的加权特征向量；共享相位的对照用于隔离特征分解器差异。
- 固定初始中心后的标签、中心、平方距离、类内距离总和及轮廓系数。
- 占比、持续时间、排序、个体转移及组平均。

必需检查失败时返回非零退出码。原始特征向量的共享相位和端到端对照均为必需检查；滤波差异另行报告。报告与具体运行版本和输入绑定；通过检查不等于 MATLAB 等价或科学结论有效。

## 3. 已记录的验证

Windows 便携环境：Octave 11.3.0、signal 1.4.7、statistics 1.8.3。误差及 Python 依赖版本见 `report_octave_11.3.0.json`。

这份历史记录使用了现已撤销的额外定号规则，不能作为当前实现的验证结果。当前版本保留求解器返回的符号，并将原始特征向量差异设为必需检查。此验证没有运行完整原始脚本或真实数据分析。

## 4. 持续集成（CI）

GitHub Actions 在主分支 push 和 PR 时运行 Python UT、Octave 对照和 MATLAB R2025b 测试。两者使用 uv 锁文件重建环境；Octave 任务安装 Linux 工具包和无窗口绘图后端，失败时也上传中间结果、日志和报告。

MATLAB 任务在已有 Python 调用测试之后，复用交接包运行独立的 `reference.m`（需要 Signal Processing Toolbox 和 Statistics and Machine Learning Toolbox），再由 Python 读取 `reference.mat` 对照生产函数。预处理直接调用 `preprocessing.py`，不在比较器中另写一份滤波实现。

MATLAB CI 使用 `--strict-signal`，将滤波信号和相位误差也设为必需检查（绝对容差 1e-8）；原始特征向量对照始终为必需检查。聚类采用共享固定初始中心，状态统计采用预设标签，不能据此推断随机聚类或论文全部结果等价。输入、参考脚本、日志和逐项误差报告作为 `matlab-reference-report` 上传。

原始特征向量不匹配会使 CI 失败；不通过修改两边符号规则来获得绿色结果。Linux 与 Windows 的 Octave 版本可能不同，因此日志会记录运行环境。

## 5. 委托别人用 MATLAB 测试

不需要本机安装 MATLAB，即可准备交接包：

```bash
uv run --locked python validation/compare.py --prepare-only --output validation_output/matlab_handoff
```

将整个目录压缩交给操作者。包内的 `HANDOFF.md` 包含中文步骤；对方只需运行 `run_reference`，回传目录及生成的 `reference.mat`、`runtime.log`。输入、脚本校验值和环境日志随包保留。同一包也可以先用 Octave 演练。

收到结果后，在无需 MATLAB 的环境中运行：

```bash
uv run --locked python validation/compare.py --compare-only --output /path/to/returned-folder
```

请在本地保留原始清单并与回传清单核对。SHA256 清单用于检测意外修改，不是防止蓄意替换的数字签名。此流程对照独立参考实现；拿到实际 MATLAB 回传结果前，不能宣称已完成 MATLAB 验证。

## 在 MATLAB 中直接测 Python

从仓库根目录执行 `uv sync --locked`，然后用 MATLAB R2025b 打开同一目录，执行：

```matlab
results = runtests('tests/matlab');
disp(table(results));
assertSuccess(results);
```

测试使用项目 `.venv` 的 Python 3.11，在内存中构造输入。每个用例直接调用生产函数并检查手算结果，无需 `input.mat`、受试者数据或额外 MATLAB 工具箱。`validation/` 的双实现数值对照是另一套检验。
