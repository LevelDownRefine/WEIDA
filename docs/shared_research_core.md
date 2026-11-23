# 共用研究接口：先在 WEIDA 内抽取

## 本次完成

| 层 | 接口 | 约定 |
| --- | --- | --- |
| 数据读取 | `research_core.io.load_scans`、`read_saved_scan` | 脑区 × 时间；保留扫描身份、条件、脑区顺序；检查患者配对 |
| 幅度预处理 | `research_core.preprocessing.preprocess_bold` | 去线性趋势、去均值、二阶 Butterworth 滤波；返回幅度序列 |
| 模块汇总 | `research_core.networks.module_summary` | 显式指定有向性；计算模块内、模块间、内外差和节点内外差 |
| 研究比较 | `research_core.statistics.compare_feature_families` | 每个指标族为扫描 × 特征；返回 ON−OFF、OFF−HC、ON−HC |

WEIDA 现有入口已调用这些实现。`run_weida.load_scans` 和 `downstream` 中的两个统计函数保留兼容导入，原有脚本及 MATLAB 委托接口无需改名。

读取支持 fmripost 带表头 CSV 与 PyLEiDA 的 `time_series/…/tseries.csv`。新元数据记录文件格式和绝对路径；旧元数据中的 `tseries.csv` 按 PyLEiDA 格式回读。下游检查帧数、非有限值和可用的脑区表头，避免把第一帧当表头丢失。无表头文件的脑区顺序依赖数据集标签文件，本身不能独立验证顺序。

## 算法边界

```text
读取 BOLD → 选定幅度预处理
                ├─ WEiDA：相位 → 加权特征向量 → 聚类 → 状态指标
                └─ PhiID：官方估计器 → 各分量连接矩阵 → 定义研究指标（待接入）
                                   ↓
                    扫描 × 特征 → 共用研究比较
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

保留现有 WEIDA 行为：ON/OFF 双侧配对 t 检验；与 HC 比较采用双侧等方差独立样本 t 检验；逐特征排除缺失值。每个比较、每个指标族分别对有限 p 值做 BH；不能检验的列保留 NaN。该策略与 PhiID 现有将无效槽位按 1 纳入校正的策略有区别，跨仓库接入时需显式统一，不能声称当前已经完全共用。

## 验证与后续

- 新增幅度保持和频带检查、有向边手算及节点重排检查、配对与校正族检查。
- 临时合成两种数据目录，实际运行 `run_features.py`，核对两种格式的特征和三组统计结果。
- 重构前固定随机输入快照与重构后逐项比较：加权特征向量、扫描长度、功能分离度、配对及独立统计结果完全相同。这只验证本次重构，不等于 MATLAB 等价验证。
- 新测试由现有 CI 的 unittest 自动发现，无需另外添加工作流。

尚未完成：把这个核心打包供 PhiID 直接依赖、迁移 PhiID 调用、共用作图接口、临床关联、基因及空间零模型。本次未重跑真实队列，也未产生新的显著性结论。
