# WEiDA 脑状态分析：Python 迁移版

本项目按原始 [StateSpace_modelling.m](https://github.com/chenfei-ye/WEIDA-DBS/blob/main/StateSpace_modelling.m) 的四个编号阶段组织代码：

| 阶段 | Python 函数 | 作用 |
| --- | --- | --- |
| 1 | `compute_weighted_eigenvectors` | 去趋势、去均值、0.04–0.07 Hz 滤波、Hilbert 相位及加权特征向量 |
| 2 | `cluster_brain_states` | 四状态 K-means，300 次初始化、最多 500 次迭代，计算距离、类内距离总和及轮廓系数 |
| 3 | `analyze_state_occupancy` | 计算状态占比和持续时间，排除扫描首尾不完整的状态段 |
| 4 | `calculate_state_dynamics` | 按总体出现频率排序状态，计算每次扫描及各组的转移概率 |

## 环境与运行

使用 uv 管理独立环境。Python 固定为 3.11.14，依赖声明在 `pyproject.toml`，具体版本锁定在 `uv.lock`，本地环境位于 `.venv`。

在仓库根目录执行：

```bash
uv sync --locked
uv run --locked python run_weida.py /path/to/PD_DBS_zcc/derivatives/fmripost output --exclude-stn 18 --tr 2.1
uv run --locked python -m unittest discover -s tests -t . -v
```

当前数据包含 15 个 HC 和 26 对 STN。STNon18 有 53 个全零脑区，显式指定 `--exclude-stn 18` 后使用 15 个 HC、25 个 OFF 和 25 个 ON 扫描。入口检查脑区列名、重复扫描和 ON/OFF 配对，只选择 HC 与 STN。排除项会写入元数据，输入文件不会被修改。正式分析前仍需确认实际 TR 和研究纳入标准。

CI 使用同一 Python 版本和 `uv sync --locked`。如需升级依赖，执行 `uv lock --upgrade`，重跑单元测试和 Octave 对照后再提交锁文件。

## 数值约定与限制

当前完成的是公式层面的迁移，尚未证明与 MATLAB 数值等价。

原脚本先把多个特征向量加权相加，再统一调整总向量符号。单个特征向量的符号任意，替换特征分解器可能改变结果。本版先令每个特征向量绝对值最大分量为正，再应用原脚本的总向量符号规则。这是明确的方法改动，不是已证实等价的移植；近简并特征空间仍存在方向歧义。

滤波显式使用延拓长度 `3*(max(len(a),len(b))-1)`，但尚未验证所有 MATLAB 版本的边界行为。不同库的 K-means 初始化也不同，比较结果应先匹配状态，不能直接比较状态编号。轮廓系数使用平方欧氏距离精确计算，大数据时开销较大。

## 输出文件

- `clustering.npz`：排序前的标签、聚类中心、平方距离、类内距离总和及轮廓系数。
- `results.npz`：排序后的中心、占比、持续时间、从 0 开始的状态序列、扫描转移矩阵、有效源状态标记及组平均转移矩阵。
- `metadata.json`：扫描顺序、受试者标识、脑区名称和参数；`lengths` 保存各扫描长度，便于恢复边界。
- `Vemp.csv`：脑区 × 状态；`State_ts.csv`：从 1 开始的状态序列。
- `Ptrans1.csv` 等：按元数据中的扫描顺序导出的个体转移矩阵。
- `P1emp.csv`、`LT1emp.csv`、`PTR1emp.csv` 等：编号 1/2/3 对应 HC/OFF/ON 的占比、持续时间及组转移矩阵。

占比和持续时间 CSV 只保留对应组的行，与原 MATLAB 使用全体扫描行并填充 NaN 的布局不同。单次扫描没有某源状态时，其转移行置零；组平均排除该扫描在此源状态的贡献，整组无有效贡献时为 NaN。

## 测试与委托验证

详见 [测试说明](validation/README.md) 和 [MATLAB 委托测试说明](validation/HANDOFF.md)。

现有 30 个单元测试和独立 Octave 对照均已接入 CI。便携 Octave 11.3.0 的记录中，17 项比较有 16 项通过，原符号加权特征向量不一致。这项差异单独报告，不代表已经复现论文，也不代表 MATLAB 等价验证通过。

完整真实数据运行及 MATLAB 基准对照仍未完成。委托 MATLAB 操作者无需安装 Python 或 uv；只有准备测试包、接收并比较结果的一方使用 uv。

原研究代码及论文来源见上方链接；本仓库不据此宣称获得上游原始源码的再分发授权。

## 结果绘图

```bash
uv run --locked python plot_results.py /path/to/output
```

默认在 `figures/` 中输出中文 PNG 和 PDF：状态概览、组转移概率、78 脑区状态模式和全部扫描状态序列。误差条为扫描间标准误，不代表显著性检验。中文显示需要微软雅黑、黑体或 Noto Sans CJK SC 字体。

## 下游分析模块

`downstream.py` 提供以下接口，对应上游预测、MSanalysis、calcu_FS 与组间比较的主要计算：

- `predict_bold` / `predict_states`：用已有中心分配新扫描，不重新训练或按新队列频率排序。训练与预测必须使用同一脑区顺序和符号约定。无完整段持续时间返回 NaN（上游预测函数返回 0）。
- `state_metastability`：按上游公式计算 `sum((特征-固定中心)^2)/(脑区数*(状态帧数-1))`，少于两帧返回 NaN。该指标是特征空间离散度，不是 Kuramoto 相位同步序参量的标准差。
- `functional_segregation`：按指定状态中心的正/非正权重分成两模块，用 Fisher z 相关计算模块内均值减模块间均值，以及节点指标。两个模块均需至少两个脑区；相关系数裁剪到 ±(1-1e-7)，返回裁剪边数。
- `paired_comparison`：按受试者 ID 检验 ON−OFF，拒绝重复或不完整配对。每个指标按成对有限值筛选，返回有效人数；零差值方差时不作退化 t 推断。
- `independent_comparison`：HC 与 OFF/ON 的独立样本 t 检验，默认等方差以对应原脚本，可通过接口显式选择 Welch。

```bash
uv run --locked python run_features.py /path/to/model/output --output /path/to/new/downstream_output --module-state 2
```

输入元数据中的原始 CSV 路径必须可访问。输出目录必须尚不存在，避免覆盖已有分析。模块状态编号从 1 开始；我们的状态 2 未必对应论文状态 2，因此不能直接套用论文解释。

输出包括 `features.npz`、`paired_statistics.json`、`independent_statistics.json` 和元数据。BH-FDR 分别作用于每一比较、每一指标族：四个状态、16 条转移或78个节点分别校正，不能解释为对所有分析统一控制。没有有效源状态的转移行先标为缺失再统计。JSON 中不可估计的数值目前使用 Python 的 NaN 表示。

同一数据用于拟合状态再做组间分析，当前结果仅作探索；训练集中心回代标签一致率用于实现检查，不是外部预测准确率。TI 队列的正式预测及 MATLAB 对照尚未执行。现有 Octave 对照仍覆盖原四阶段；新增模块目前依靠单元测试及真实数据试运行验证。

## Shared amplitude preprocessing package

`common/preprocessing` is the installable `weida-preprocessing` package (NumPy/SciPy only).
WEIDA's `preprocessing.py` delegates to it; downstream projects should pin a full Git commit
that has passed this repository's Python and original MATLAB preprocessing CI.

```python
from weida_preprocessing import preprocess_scan
filtered = preprocess_scan(bold, remove_trend=True, demean=True, band=(0.04, 0.07), tr=2.1)
```

Transforms are opt-in; the default returns a copy. No phase extraction or frame trimming.
[MATLAB source attribution and coverage](tests/matlab/PREPROCESSING.md).

## 共享后处理

`common/statistics.py` 提供 `weida-common`，直接封装原作者 `group_comparison.py` 的三组 t 检验和 statsmodels BH。公共包直接调用 `statistics.py` 中的正式统计函数，并与原脚本生成的独立离线基准对照，覆盖正常、NaN 和退化输入；不再添加置信区间、自动排除缺失值或替换退化结果。工程封装理由见 [函数与出处](common/README.md)。

```sh
uv run --locked python run_postprocessing.py --input /path/to/model_output --output /path/to/shared_report
```

可添加 `--features /path/to/run_features_output`。输出 CSV 与方法元数据；未对齐的通用绘图暂不提供。旧 `run_features.py` 不变，与本入口的缺失值及退化结果策略不同。
