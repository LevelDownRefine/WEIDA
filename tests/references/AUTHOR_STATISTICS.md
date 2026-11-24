# 作者统计代码整理

来源：[group_comparison.py](https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py)。仅将 DBS 段拆为函数，保留缺失值、退化检验、总体标准差及逐对比 BH 的原行为。

| 函数 | 作用 | 原位置 |
| --- | --- | --- |
| `describe_groups` | 各组均值、标准差 | 82–87 行 |
| `compare_groups` | 三组 t 检验 | 89–91 行 |
| `get_fdr` | 单个对比的 BH | 6–9 行 |
| `group_comparison` | 串联检验和校正 | 115–117 行 |

输入：已配对的 HC、OFF、ON，均为受试者 × 指标。组顺序 HC/OFF/ON，对比顺序 OFF−HC/ON−OFF/ON−HC。此层不做 ID 配对、绘图、文件读取或额外推断。

`statistics_fixture.json` 由固定版本原脚本直接执行生成，覆盖正常、NaN、常数列和恒定变化；不是用整理后的函数生成。基准保留原脚本的 16 HC、25 对、4 指标规模；函数支持可变指标数。

重建：`uv run --locked python tests/references/generate_statistics_fixture.py`（需要 gh 和网络）。日常测试仅读取基准：`uv run --locked python -m unittest tests.test_author_statistics`。

正式实现位于 `common/postprocessing/src/weida_postprocessing/author_statistics.py`；安装包直接提供这些入口。此目录仅保留原脚本基准及生成工具，不再复制统计实现。
