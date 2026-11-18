# 单元测试、Octave 对照与 MATLAB 委托测试

所有命令在仓库根目录运行。测试使用合成数据，不需要上传真实受试者数据。

## 1. 单元测试（UT）

```bash
uv run --locked python -m unittest discover -s tests -t . -v
```

现有 22 个用例覆盖：完整状态段的持续时间、缺失状态和全程单一状态、扫描边界隔离、按扫描平均的组转移概率、状态排序一致性、非法长度与标签、全零及非有限输入、不同扫描长度、交接包生成和校验。

预期值来自可手算样例；测试不依赖 MATLAB 或 Octave。

## 2. 当前的 Octave 对照

```bash
uv run --locked python validation/compare.py --octave C:/path/to/octave-cli.exe --output validation_output
```

比较器生成固定种子的 78 脑区、210 时间点数据，调用 `reference.m`，输出中间 `.mat` 文件、日志和 `report.json`。

参考脚本按原公式独立实现，并非未经修改的整份上游脚本：数据读取改为固定测试输入；分别导出原符号行为和统一符号后的结果；聚类使用分离的模拟簇与共享初始中心，而非随机重启 300 次；状态统计使用预设标签，涵盖缺失状态、全程单一状态及首尾不完整段。

检查内容：

- 滤波系数、滤波信号和相位。
- 原符号与统一符号的加权特征向量；共享相位的对照用于隔离特征分解器差异。
- 固定初始中心后的标签、中心、平方距离、类内距离总和及轮廓系数。
- 占比、持续时间、排序、个体转移及组平均。

必需检查失败时返回非零退出码。滤波和原符号端到端差异另行报告，不强制相等。报告与具体运行版本和输入绑定；通过检查不等于 MATLAB 等价或科学结论有效。

## 3. 已记录的验证

Windows 便携环境：Octave 11.3.0、signal 1.4.7、statistics 1.8.3。误差及 Python 依赖版本见 `report_octave_11.3.0.json`。

17 项比较中 16 项通过，原符号特征向量仍存在不一致，说明统一符号是方法变更。此验证没有运行完整原始脚本或真实数据分析。

## 4. 持续集成（CI）

GitHub Actions 在 push 和 PR 时分别运行 Python UT 与 Octave 对照。两者使用 uv 锁文件重建环境；Octave 任务安装 Linux 工具包和无窗口绘图后端，失败时也上传中间结果、日志和报告。

CI 绿色仅表示必需检查通过，不能掩盖报告中的原符号差异。Linux 与 Windows 的 Octave 版本可能不同，因此日志会记录运行环境。

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
