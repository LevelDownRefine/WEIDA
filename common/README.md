# 共用前后处理

## 已完成什么

- 将共用代码集中到根目录 `common/`，预处理和后处理分别打包、独立安装。
- 预处理实现从算法中独立，WEIDA 调用同一份实现，避免重复滤波。
- 将已有通用后处理抽成共享包，保留指标单位、定义及预先声明的校正族。
- 增加 WEIDA 已保存结果的适配器和 `run_postprocessing.py` 入口，无需重新拟合状态。
- 同步 uv workspace、锁文件及安装路径；调用方使用已通过 CI 的固定提交。

## 目录与职责

| 目录 | 包名 / 导入名 | 职责 |
| --- | --- | --- |
| `preprocessing/` | `weida-preprocessing` / `weida_preprocessing` | 脑区 BOLD 幅值预处理 |
| `postprocessing/` | `weida-postprocessing` / `weida_postprocessing` | 指标配对、组间统计、校正、通用作图及保存 |

两个包都不包含 WEiDA 相位/特征分解/聚类或 PhiID 分解算法。状态占比、持续时间、转移概率等算法专属指标不在这里计算；结果适配器只读取已经计算好的指标。

## 预处理

入口是 `preprocess_scan`，输入和输出均为 **脑区 × 时间点**：

```python
from weida_preprocessing import preprocess_scan

filtered = preprocess_scan(
    bold, remove_trend=True, demean=True, band=(0.04, 0.07), tr=2.1,
)
```

按时间轴对每个脑区独立执行：

1. `detrend`：去线性趋势。
2. 减去时间均值。
3. 二阶 Butterworth 带通及 `filtfilt` 双向滤波。

默认调用 `preprocess_scan(bold)` 仅返回副本，变换需要显式启用。输入不被修改，输出不裁剪首尾、不提取相位。

## 后处理

输入是 `FeatureTable`：扫描 × 指标的数值，配合受试者/条件、指标名、单位、方法说明和校正族。

| 入口 | 功能 |
| --- | --- |
| `group_indices` | 按受试者 ID 对齐 ON/OFF，检查重复记录和配对完整性 |
| `compare_features` | ON−OFF 配对双侧 t 检验，OFF/ON−HC 等方差独立双侧 t 检验；返回均差、有效人数、未校正 95% 置信区间和 p/q |
| `adjust_bh` | 每个对比 × 指标族分别做 BH 校正 |
| `plot_features` | 按指标族绘制分布/散点图与患者配对线，输出 PNG/PDF |
| `write_table` | 保存 CSV 表格 |
| `save_results` | 统一保存统计表、个体指标表、元数据和图件，拒绝覆盖已有目录 |
| `adapters.load_weida` | 将已有 WEIDA 输出转为指标表，保留状态编号及各指标单位 |

```python
from weida_postprocessing import compare_features, save_results

statistics = compare_features(table)
save_results(table, 'results/shared_report')
```

WEIDA 命令行入口：

```sh
uv run --locked python run_postprocessing.py --input /path/to/model_output --output results/shared_report
```

可添加 `--features /path/to/run_features_output`。输出包括 `statistics.csv`、`subject_features.csv`、`metadata.json` 和图件。

配对比较成对排除缺失值，独立比较按组/指标排除。不可估计的检验结果为 NaN，但保留预先定义的 BH 族规模。旧 `downstream.py` 只对有限 p 值校正，两者在缺失/退化指标上可能不同；旧入口的行为没有改动。

## 已验证与尚未覆盖

- 预处理：真实 MATLAB R2025b 直接调用封装的上游参考函数，与安装后的 Python 包对照；覆盖 78 脑区、4 种长度和不同 TR，包含首尾时间点。参考函数的中文注释注明固定版本与出处，见 [预处理对照](../tests/matlab/PREPROCESSING.md)。
- 后处理：Python 测试覆盖统计数值、ID 重排、缺失值、BH 族规模及 CSV/PNG/PDF 输出，见 [后处理测试](../tests/test_postprocessing.py)。
- 上游统计对照另由 [PR #13](https://github.com/LevelDownRefine/WEIDA/pull/13) 维护；它验证现有 `downstream.py`，尚不能代替本共享包的上游等价验证。
- 不宣称真实队列、论文结果、绘图外观或全部后处理已与原作者完全等价。
- PhiID 的预处理依赖已接入；共享后处理的消费端接入仍待推进。

## 安装与维护

仓库内执行 `uv sync --locked`，再运行 `uv run --locked python -m unittest discover -s tests -t .`。

外部调用方固定 WEIDA 完整提交 SHA，并选择 `common/preprocessing` 或 `common/postprocessing` 子目录安装；只安装需要的包。更新依赖前核验对应提交的 CI：预处理需 Python 与 MATLAB 对照，后处理需其相应测试。测试统一维护在根目录 `tests/`，不在调用方复制实现。
