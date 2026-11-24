# 共用前后处理

## 预处理（`weida_preprocessing`）

| 函数 | 做什么 | 出处 |
| --- | --- | --- |
| `preprocess_scan` | 按需去线性趋势、去均值、带通及双向滤波，保留脑区与时间点 | 根据原作者 [StateSpace_modelling.m 第 70–88 行](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L70-L88) 的幅值处理段封装 |

## 后处理（`weida_postprocessing`）

| 函数 / 类型 | 做什么 | 出处 |
| --- | --- | --- |
| `FeatureTable` | 保存扫描 × 指标及其名称、单位、校正族 | 自编，抽自 PhiID 的 `research.py` |
| `group_indices` | 按受试者 ID 配对 ON/OFF，检查重复与缺失 | 自编，抽自 PhiID 的 `research.py` |
| `compare_features` | 三组 t 检验、均差、95% 置信区间及 p/q | 自编；检验方式参考原作者 [group_comparison.py](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)，置信区间为新增 |
| `adjust_bh` | 按对比和指标族做 BH 校正 | 自编；参考同上，额外保留不可估计项的族位置 |
| `plot_features` | 分布、散点和患者配对图，输出 PNG/PDF | 自编，抽自 PhiID 的 `research.py` |
| `write_table` | 保存 CSV 表 | 自编，抽自 PhiID 的 `research.py` |
| `save_results` | 保存统计表、个体指标、元数据及图件 | 自编，抽自 PhiID 的 `scripts/analyze_results.py` |
| `adapters.load_weida` | 将已有 WEIDA 结果转为 `FeatureTable` | 自编，抽自 PhiID 的 `result_adapters.py` |

后处理是自编公共实现；参考原作者的方法不等于已验证全部输出与其完全一致。
