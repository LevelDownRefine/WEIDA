# 统一研究入口：WEiDA / PhiID

## 总接口

```python
from research_core.study import run_study

for method in ("weida", "phiid"):
    result = run_study(
        "/path/to/input", f"output/{method}",
        method=method, preprocessing="weida_amplitude", tr=2.1,
        exclude_stn=["18"],
    )
```

两种算法走同一条读取、预处理、统计、保存和作图路径。已有内存数据可调用 `analyze(scans, metadata, regions, method=..., preprocessing=...)`。算法专属参数通过 `options` 传入：WEiDA 为 `k/seed/replicates`，PhiID 为 `lag`（帧数）。不适用的参数报错，不会静默忽略。

`preprocessing` 必须显式传入。`loaded` 保留输入幅度；`weida_amplitude` 采用现有去趋势、去均值及带通处理。WEiDA 后端只做后续相位和特征向量步骤，不再次滤波；PhiID 后端调用官方 `phyid.calculate.calc_PhiID`，依赖提交固定在 `uv.lock`，与 PhiID 仓库采用同一官方版本。不依赖本机相邻仓库或 Python 路径注入。

## 统一输出

`StudyResult` 包含算法结果、扫描元数据、脑区顺序、预处理记录和三组比较。算法返回 `AlgorithmResult`：指标族、列名、单位、定义、参数和原始结果。两者都通过同一个统计函数与图形函数。

| 文件 | 内容 |
| --- | --- |
| `features.npz` | 各指标族的扫描 × 特征数组 |
| `artifacts.npz` | 算法原始结果；WEiDA 状态/中心/向量，PhiID 完整连接矩阵及上三角存储 |
| `metadata.json` | 方法、输入数值指纹、扫描顺序、脑区标签、参数、指标含义、单位和统计约定 |
| `statistics.json` | ON−OFF、OFF−HC、ON−HC，各指标族分别校正；无效统计量保存为 null |
| `figures/` | 两种算法共用的组分布图与患者 ON/OFF 配对变化图，PNG/PDF |

WEiDA 默认返回占比、持续时间、转移概率和状态离散度。PhiID 默认返回 16 分量全连接均值，纳入每条边的两个端点方向；同时保留扫描 × 脑区 × 脑区 × 16 的矩阵，对角线为 NaN。**这套默认检验是全连接均值检验，不是逐边、NBS 或整图 RMS 检验。** 后续的节点/子网络指标需给出明确的汇总定义后接入。两种算法的指标含义并未合并。

原 `run_weida.py` / `run_features.py` 仍使用原输出文件名，兼容已有使用者；新入口输出遵循以上统一格式。

## 本次完成

| 层 | 接口 | 约定 |
| --- | --- | --- |
| 数据读取 | `research_core.io.load_scans`、`read_saved_scan` | 脑区 × 时间；保留扫描身份、条件、脑区顺序；检查患者配对 |
| 幅度预处理 | `research_core.preprocessing.preprocess_bold` | 去线性趋势、去均值、二阶 Butterworth 滤波；返回幅度序列 |
| 模块汇总 | `research_core.networks.module_summary` | 显式指定有向性；计算模块内、模块间、内外差和节点内外差 |
| 研究比较 | `research_core.statistics.compare_feature_families` | 每个指标族为扫描 × 特征；返回 ON−OFF、OFF−HC、ON−HC |
| 总入口 | `research_core.study.run_study` / `analyze` | 选择后端，共用完整研究流程 |
| 作图 | `research_core.plotting.plot_study` | 从统一结果读取数值、名称和单位 |

WEIDA 现有入口已调用这些实现。`run_weida.load_scans` 和 `downstream` 中的两个统计函数保留兼容导入，原有脚本及 MATLAB 委托接口无需改名。

读取支持 fmripost 带表头 CSV 与 PyLEiDA 的 `time_series/…/tseries.csv`。新元数据记录文件格式和绝对路径；旧元数据中的 `tseries.csv` 按 PyLEiDA 格式回读。下游检查帧数、非有限值和可用的脑区表头，避免把第一帧当表头丢失。无表头文件的脑区顺序依赖数据集标签文件，本身不能独立验证顺序。

## 算法边界

```text
读取 BOLD → 选定幅度预处理
                ├─ WEiDA：相位 → 加权特征向量 → 聚类 → 状态指标
                └─ PhiID：官方估计器 → 各分量连接矩阵 → 全连接均值
                                   ↓
                    扫描 × 特征 → 共用研究比较 → 保存 / 作图
```

- Hilbert 相位、特征向量、聚类、状态占比及离散度保留在 WEiDA 算法层。
- WEiDA 功能分离度仍由 `functional_segregation` 包装：中心正/非正定义模块，Pearson 相关转 Fisher z，再调用模块汇总。
- 该功能分离度按原入口使用读入的幅度序列；本次没有偷偷增加一次带通滤波。
- PhiID 矩阵不可套用 Fisher z。调用共用模块汇总前，需单独确定分量及模块定义。
- 有向矩阵必须使用两个方向的全部非对角元素；节点值表示出向连接。无向矩阵要求对称，每条边计一次。不能默认为 PhiID 单个分量取上三角，也不能随意混合 X/Y 分量。

## 统计接口

```python
from research_core.statistics import compare_feature_families

comparisons = compare_feature_families(
    features={"occupancy": occupancy, "lifetime": lifetime},
    subject_ids=participant_ids,
    conditions=conditions,
)
```

输入每一行对应一次扫描；患者 ON/OFF 使用同一受试者 ID。HC 使用独立 ID。特征列可以是状态、节点或预先定义的连接指标，但统计单位始终是受试者。

统一入口对两个后端都采用现有 WEIDA 行为：ON/OFF 双侧配对 t 检验；与 HC 比较采用双侧等方差独立样本 t 检验；逐特征排除缺失值。每个比较、每个指标族分别对有限 p 值做 BH；不能检验的列保留 NaN。PhiID 仓库旧独立入口仍将无效槽位按 1 纳入校正，因此旧入口的 q 值不能在有无效列时直接混用。本 PR 内的两个后端已经使用同一个统计实现。

## 验证与后续

- 新增幅度保持和频带检查、有向边手算及节点重排检查、配对与校正族检查。
- 临时合成两种数据目录，实际运行 `run_features.py`，核对两种格式的特征和三组统计结果。
- 重构前固定随机输入快照与重构后逐项比较：加权特征向量、扫描长度、功能分离度、配对及独立统计结果完全相同。这只验证本次重构，不等于 MATLAB 等价验证。
- 新测试由现有 CI 的 unittest 自动发现，无需另外添加工作流。
- 两个真实算法均经同一 Python API 和 CLI 运行测试，无 mock 估计器；核对官方 PhiID 数值及 bit 单位、脑区重排和 X/Y 交换、WEiDA 旧接口一致性与只滤波一次、共用输出及 PNG/PDF。

尚未完成：迁移 PhiID 仓库的旧独立入口以消除跨仓库重复适配、统一总入口中的模块指标配置、临床关联、基因及空间零模型。本次未重跑真实队列，也未产生新的显著性结论。
