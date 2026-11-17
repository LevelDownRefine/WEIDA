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
uv run --locked python -m unittest discover -v
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

现有 9 个单元测试和独立 Octave 对照均已接入 CI。便携 Octave 11.3.0 的记录中，17 项比较有 16 项通过，原符号加权特征向量不一致。这项差异单独报告，不代表已经复现论文，也不代表 MATLAB 等价验证通过。

完整真实数据运行及 MATLAB 基准对照仍未完成。委托 MATLAB 操作者无需安装 Python 或 uv；只有准备测试包、接收并比较结果的一方使用 uv。

原研究代码及论文来源见上方链接；本仓库不据此宣称获得上游原始源码的再分发授权。
