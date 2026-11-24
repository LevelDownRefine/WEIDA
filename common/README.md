# 共用前后处理

## 预处理（`weida_preprocessing`）

| 函数 | 做什么 | 出处 |
| --- | --- | --- |
| `preprocess_scan` | 预处理入口 | [原流程](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L88) |
| `scipy.signal.detrend` | 去线性趋势 | [detrend](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L87) |
| `ndarray.mean` | 求均值后减去 | [demean](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/functions/demean.m) |
| `scipy.signal.butter` | 构造带通滤波器 | [butter](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L75) |
| `scipy.signal.filtfilt` | 双向滤波 | [filtfilt](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L88) |

## 后处理（`weida_postprocessing`）

| 函数 / 类型 | 做什么 | 出处 / 封装理由 |
| --- | --- | --- |
| `compare_features` | 三组 t 检验、均值/标准差、均差及 p/q | [作者检验](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)；均差仅汇总效应方向 |
| `adjust_bh` | statsmodels BH 校正 | [get_fdr](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py) |
| `FeatureTable` | 指标、单位和校正族 | 工程封装：避免混淆指标与检验族 |
| `group_indices` | ID 配对与完整性检查 | 工程封装：避免按文件顺序错配患者 |
| `write_table` | 保存 CSV | 工程封装：统一输出字段 |
| `save_results` | 保存统计、个体指标和元数据 | 工程封装：记录方法、避免覆盖 |
| `adapters.load_weida` | 读取已有 WEIDA 结果 | 工程封装：接入已有输出，不重算指标 |

统计不删除 NaN、不替换退化结果、不增加置信区间；每个对比按预先声明的指标族校正。直接对照公共包的测试见 [test_postprocessing.py](../tests/test_postprocessing.py)。作图暂不纳入，待对齐作者实现。
