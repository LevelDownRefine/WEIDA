# 共用前后处理

## 预处理（`weida_preprocessing`）

| 函数 | 做什么 | 出处 |
| --- | --- | --- |
| `preprocess_scan` | 串联幅值预处理，保留脑区与时间点 | 我们封装的 [Python 入口](preprocessing/src/weida_preprocessing/__init__.py) |
| `scipy.signal.detrend` | 对每个脑区去线性趋势 | 对应原版 MATLAB `detrend`：[第 87 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L87) |
| `numpy.ndarray.mean` | 计算时间均值，再从信号中减去 | 对应作者的 [`functions/demean.m`](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m)，由主脚本[第 87 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L87)调用 |
| `scipy.signal.butter` | 构造二阶 Butterworth 带通滤波器 | 对应 MATLAB Signal Processing Toolbox 的 `butter`：[第 70–75 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L75) |
| `scipy.signal.filtfilt` | 沿时间轴双向滤波 | 对应 MATLAB Signal Processing Toolbox 的 `filtfilt`：[第 88 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L88) |
| `original_preprocess`（测试） | 生成独立 MATLAB 参考结果 | 我们将原作者[第 70–88 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L88)中的幅值处理段封装为[参考函数](../tests/matlab/original_preprocess.m)，不包含 Hilbert 相位提取 |
| `demean`（测试） | MATLAB 参考代码中去均值 | 直接取自作者的 [`functions/demean.m`](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m)，仅翻译注释并补充出处 |

## 后处理（`weida_postprocessing`）

| 函数 / 类型 | 做什么 | 出处 |
| --- | --- | --- |
| `FeatureTable` | 保存扫描 × 指标及其名称、单位、校正族 | 自编，抽自 PhiID 的 `research.py` |
| `group_indices` | 按受试者 ID 配对 ON/OFF，检查重复与缺失 | 自编，抽自 PhiID 的 `research.py` |
| `compare_features` | 三组 t 检验、均差、95% 置信区间及 p/q | 自编；检验方式参考原作者 [group_comparison.py](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)，置信区间为新增 |
| `adjust_bh` | 按对比和指标族做 BH 校正 | 自编，调用 SciPy；参考同上，额外保留不可估计项的族位置 |
| `plot_features` | 分布、散点和患者配对图，输出 PNG/PDF | 自编，抽自 PhiID 的 `research.py` |
| `write_table` | 保存 CSV 表 | 自编，抽自 PhiID 的 `research.py` |
| `save_results` | 保存统计表、个体指标、元数据及图件 | 自编，抽自 PhiID 的 `scripts/analyze_results.py` |
| `adapters.load_weida` | 将已有 WEIDA 结果转为 `FeatureTable` | 自编，抽自 PhiID 的 `result_adapters.py` |

后处理是自编公共实现；参考原作者的方法不等于已验证全部输出与其完全一致。
