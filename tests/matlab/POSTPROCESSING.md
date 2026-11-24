# 后处理上游对照

参考版本统一为 WEIDA-DBS `2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c`。参考函数中文注释包含原文件链接；只封装计算段和输入输出，不保存完整脚本副本。生产算法不改动，本 PR 不包含共享后处理包。

## MATLAB 对照

| 参考函数 | 原作者出处 | Python 对象与覆盖范围 |
| --- | --- | --- |
| `original_state_occupancy` | `StateSpace_modelling.m` 第 3 节 | `analyze_state_occupancy`：占比、以秒计的完整段持续时间、首尾不完整段、缺席状态及单帧扫描。 |
| `original_predict_states` | `functions/WEiDA_fix_cluster.m` 第 2 节 | `predict_states` 的标签：固定特征和中心，保留中心顺序，平距选第一个。 |
| `original_transitions` | 同文件末尾转移概率段 | `calculate_state_dynamics` 的每扫描转移矩阵及频率排序；包括缺席源状态和扫描边界。 |
| `original_metastability` | `MSanalysis.m` DBS 段 | `state_metastability`：78 脑区、4 状态、两次不同长度扫描的离散度。 |

固定输入特征/标签，不重新运行 Hilbert、eigs 或 kmeans。有限数值绝对容差 `1e-12`，标签和缺失掩码精确比较。MATLAB 测试在现有 R2025b CI 运行，不依赖 Octave。

## 作者 Python 对照

`tests/references/original_postprocessing.py` 保存上游计算函数/语句的封装，`tests/test_upstream_postprocessing.py` 执行：

- `group_comparison.py`：HC 16 人、OFF/ON 25 对、4 指标的配对/等方差独立 t 检验及 statsmodels BH，与现有 downstream 的 t/p/q 对照；同时检查 ID 重排。
- `calcu_FS.py`：模块内外相关经 Fisher-z 后的差值、逐节点分离度，覆盖 4/7/78 脑区和不等大小模块。
- 参考函数依赖 pandas、statsmodels，仅加入开发测试依赖；数值使用 `atol=1e-12, rtol=0`。

## 明确的边界

- 持续时间使用主脚本的 NaN 规则；`WEiDA_fix_cluster.m` 无完整段时返回 0，不能声称其全部输出与当前 Python 等价。
- 原离散度公式在单帧状态可能产生 Inf、缺席状态返回 0；当前 Python 对不足两帧返回 NaN。单独测试这项既有差异，不改容差或隐藏结果。
- 原组间脚本不处理缺失值；当前 Python 会成对/按指标排除。单独记录缺失输入的不同结果；正常有限输入才宣称数值对齐。
- 完美相关的 Fisher-z 原公式未裁剪；当前 Python 显式裁剪。功能分离度等价测试使用无需裁剪的输入。
- 网络重叠暂无主分支 Python 实现，临床关联尚不在该仓库；本测试 PR 不新增这些算法，也不宣称验证绘图外观或论文队列结果。

运行：`uv run --locked python -m unittest discover -s tests -t .`；MATLAB 使用现有 `runtests('tests/matlab')`。
