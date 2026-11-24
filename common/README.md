# 正式统计入口

`from common.statistics import describe_groups, compare_groups, get_fdr, group_comparison`

| 函数 | 作用 |
| --- | --- |
| `describe_groups` | 各组均值和总体标准差 |
| `compare_groups` | OFF−HC、ON−OFF、ON−HC 的 t/p |
| `get_fdr` | 单个对比的 BH 校正 |
| `group_comparison` | 串联检验与校正，返回 t/p/q |

输入为已配对的 HC、OFF、ON 受试者 × 指标矩阵。实现位于 `statistics.py`，保留[作者原行为](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)。测试直接调用本包，使用原脚本生成的独立基准。

WEiDA 和 PhiID 均可传入受试者 × 指标矩阵，共用这套统计；ON/OFF 须按患者配对，HC 为独立组。

## PhiID 共用读写

| 函数 | 用途 | 出处 |
| --- | --- | --- |
| `read_bold_csv` | 去索引列，转为脑区 × 时间 | [importdata](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L25-L39) |
| `load_scans` | 加载扫描、标签、组别 | 原读取流程的适配；支持现有 fmripost/PyLEiDA 布局，检查标签和数据 |
| `group_indices` | 按 ID 对齐 ON/OFF | 工程补充：避免固定人数和文件顺序导致错配 |
| `write_table` | 保存带表头的结果 CSV | 工程补充：PhiID 需保存分量名、组别及数值；MATLAB readtable 验证互操作 |

实现均在 `common/io.py`。读取对照作者 MATLAB 原语句；表格保存保留 Python 双精度，不复制 csvwrite 的旧数值格式。NPZ、JSON 及算法专属结果组织仍由调用方维护。
