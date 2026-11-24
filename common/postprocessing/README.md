# 正式统计入口

`from weida_postprocessing import describe_groups, compare_groups, get_fdr, group_comparison`

| 函数 | 作用 |
| --- | --- |
| `describe_groups` | 各组均值和总体标准差 |
| `compare_groups` | OFF−HC、ON−OFF、ON−HC 的 t/p |
| `get_fdr` | 单个对比的 BH 校正 |
| `group_comparison` | 串联检验与校正，返回 t/p/q |

输入为已配对的 HC、OFF、ON 受试者 × 指标矩阵。实现位于 `src/weida_postprocessing/statistics.py`，保留[作者原行为](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)。测试直接调用本包，使用原脚本生成的独立基准。

WEiDA 和 PhiID 均可传入受试者 × 指标矩阵，共用这套统计；ON/OFF 须按患者配对，HC 为独立组。
