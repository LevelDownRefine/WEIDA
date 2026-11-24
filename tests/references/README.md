# PhiID 所需的组间统计对照

参考代码来自 WEIDA-DBS 固定提交 `2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c` 的 `group_comparison.py`，出处链接写在参考函数注释中。原作者此部分使用 Python，不需要新增 MATLAB 测试。

只保留：

- ON−OFF 配对双侧 t 检验、OFF/ON−HC 等方差独立双侧 t 检验。
- 每次对比的 BH 校正；逐项比较 t、p、q，容差 `atol=1e-12, rtol=0`。
- 合成 16 HC、25 对患者、16 个指标；检查受试者顺序打乱后配对结果不变。
- 检查 BH 对零、重复及边界 p 值的处理。
- 原脚本对 NaN 不作排除，当前 Python 成对排除；该既有差异单独记录，不宣称缺失输入完全等价。

测试对象是现有 `downstream.paired_comparison`、`independent_comparison`；尚未验证草稿共享后处理包或真实队列结果。此 PR 不新增状态占比、持续时间、转移概率、状态预测、状态离散度或功能分离度测试，也不改生产算法。

运行：`uv run --locked python -m unittest tests.test_upstream_postprocessing -v`。
