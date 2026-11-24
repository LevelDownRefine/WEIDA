# MATLAB 预处理对照

## 原版算法代码出处

固定提交：`2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c`。下列代码在 CI 中从原仓库检出，未改写。

| 代码 | 出处 | 用途 |
| --- | --- | --- |
| `fnq`、`flp`、`fhi`、`Wn`、`k`、`butter` | [StateSpace_modelling.m 第 70–76 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L76) | 根据 TR 计算奈奎斯特频率，构建 0.04–0.07 Hz、阶数参数为 2 的 Butterworth 带通滤波器。 |
| `ts=demean(detrend(BOLD(seed,:)));` | [同文件第 87 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L87) | 对每个脑区的时间序列先去线性趋势，再去均值。 |
| `signal_filt =filtfilt(bfilt,afilt,ts);` | [同文件第 88 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L88) | 双向滤波，得到进入算法前的 BOLD 幅值序列。 |
| `demean.m` | [functions/demean.m](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m) | 原作者仓库附带的去均值函数；对这里的行向量减去时间均值。 |

`detrend` 是 MATLAB 函数，`butter`、`filtfilt` 由 Signal Processing Toolbox 提供；它们不是我们实现的替代函数。原脚本紧接着的第 89 行 `angle(hilbert(...))` 属于相位提取，本测试不执行。

## 我们编写的测试代码

来源为本 PR 的 [testPreprocessing.m](testPreprocessing.m)，不是原作者算法代码。

| 函数 | 用途 |
| --- | --- |
| `testPreprocessing` | 注册 MATLAB 单元测试。 |
| `setupOnce` | 提取原版代码块，确认上游 `demean.m` 路径，连接 uv 环境中的生产 Python 入口。 |
| `testOriginalAmplitudePipeline` | 构造合成数据，执行原版预处理及 Python 预处理，逐元素比较并打印误差。 |
| `numbers` | 将 Python 数组转换为 MATLAB 数组，仅用于比较。 |
| `teardownOnce` | 恢复 MATLAB 搜索路径。 |

## 覆盖范围与运行

- 基准直接来自 [WEIDA-DBS 固定版本](https://github.com/chenfei-ye/WEIDA-DBS/tree/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c)。CI 检出该版本，不自行重写 MATLAB 参考算法。
- 从 `StateSpace_modelling.m` 原样提取滤波器设置及 `ts=demean(detrend(BOLD(seed,:))); signal_filt =filtfilt(bfilt,afilt,ts);`；`demean.m` 直接使用该仓库的 `functions/demean.m`。提取不到唯一代码块时测试报错。
- 对照安装包入口 `weida_preprocessing.preprocess_scan`：去线性趋势、去均值、二阶 Butterworth 0.04–0.07 Hz 带通及 `filtfilt`。
- 78 脑区，13/48/210/312 帧，TR 为 2.1/1.5/2.1/2.0 秒；包含所有首尾时间点，检查形状、有限值及输入不被修改。逐元素绝对容差 `1e-8`，打印每组最大误差。
- 固定种子的合成数据；不验证患者数据、fMRIPrep、Hilbert 相位、WEiDA 特征分解或 PhiID 统计结果。原脚本的加载、相位及后续算法不执行。

本地需要 MATLAB R2025b 和 Signal Processing Toolbox。在仓库根目录执行：

```sh
uv sync --locked --python 3.11
git clone https://github.com/chenfei-ye/WEIDA-DBS.git .reference/WEIDA-DBS
git -C .reference/WEIDA-DBS checkout 2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c
matlab -batch "results = runtests('tests/matlab'); disp(table(results)); assertSuccess(results)"
```

没有 MATLAB 时仍可运行现有 Python 单测；PR 中独立的 MATLAB 工作流负责真实 MATLAB 对照。
