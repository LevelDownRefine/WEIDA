# CI 数值差异调查（2026-11-23）

## 结论

- Python 的 `numpy.linalg.eigh` 替换 MATLAB 原脚本的 `eigs`，目前不能保持原始 WEiDA 特征一致；这是尚未解决的兼容性问题。
- 原脚本直接计算 `sum(lambda_i * v_i) / sum(lambda_i)`，然后只调整总向量符号。单个特征向量的相对符号变化会改变结果，不能用最终整体翻转消除。
- 不改原算法、不统一单个特征向量符号、不放宽容差、不跳过失败检查。相关 PR 保持 draft，不能宣称完整 MATLAB 等价。

## CI 实际检查

来源：[MATLAB R2025b Update 7 的运行记录](https://github.com/LevelDownRefine/WEIDA/actions/runs/45189326673)。所有输入均为合成数据。

| 检查 | 结果 |
|---|---|
| Python 单测 / Ruff | 41 项通过 / 通过 |
| MATLAB 单元与交接集成测试 | 10 项通过 |
| 严格跨语言数值对照 | 17 项中 15 项通过 |
| `original_features_shared_phase` | 失败，最大绝对差约 0.155403，容差 1e-8 |
| `original_features_end_to_end` | 失败，最大绝对差约 0.155403，容差 1e-8 |

通过项包括预处理、相位、同一组特征对的加权、固定初始中心的聚类及状态统计。特征生成失败足以阻止完整流程等价的结论。

## 隔离求解器后的证据

将 MATLAB 实际使用的 21 个 78×78 矩阵直接导出，避免在 Python 中重新计算 `cos` 带来的舍入差异：

| 比较（同一输入矩阵） | 匹配数量 | 最大绝对差 |
|---|---:|---:|
| MATLAB 重复调用 / 改全局随机种子 | 21/21 | 0 |
| MATLAB 矩阵 / 函数句柄 / 显式原起点 | 21/21 | 0 |
| NumPy `eigh` 与 MATLAB | 6/21 | 0.169306 |
| SciPy `eigsh`、使用 MATLAB 原起点（Linux CI） | 11/21 | 0.168078 |
| Python Lanczos + 三对角 QR、使用 MATLAB 原起点 | 21/21 | 2.36e-16 |
| Octave 11.3.0 与 MATLAB（Windows 本机） | 10/21 | 0.168078 |

- MATLAB 全部收敛。NumPy 特征值差与特征对残差约 1e-14；主要特征向量相关为 ±1。差异不是求解失败或特征值算错。
- MATLAB 调用轨迹显示 Krylov–Schur；其第一次数乘是输入检查，第二次才是迭代起点。
- 此版本的起点匹配 `RandStream('dsfmt19937','Seed',0,'NormalTransform','Ziggurat')`，归一化后误差约 5.55e-17。直接用 Python `randomgen.DSFMT(0)` 加 NumPy 正态采样不匹配，不能据此声称已重现初始化。
- 独立 Python Lanczos 探索仍依赖 MATLAB 导出的起点，且不是完整 Krylov–Schur 实现；仅为定位根因，不进入生产实现。
- SciPy 结果还依赖平台；上表明确记录 Linux CI，不能把 Windows 本机结果混为同一组测量。

## 独立与边界验证

[扩展验证运行](https://github.com/LevelDownRefine/WEIDA/actions/runs/45192029701)（提交 `0a884bd`）新增 100 组固定种子的独立随机相位，以及 5 个边界矩阵：恒定零相位、恒定非零相位、近恒定相位、两个相等主特征值的分块和圆周相位。

- MATLAB 的重复调用、全局种子变动、函数句柄和显式原起点均为 **126/126** 匹配，全部收敛。
- Python 探索性 Lanczos 在原 21 例和新增 100 例上匹配；新增随机组最大误差 **7.86e-16**。
- **5 个边界案例全部失败**。恒定类案例出现严重数值失稳；重根两例差异约 **0.155 / 0.160**。这是诊断候选自身的限制，不是生产 `eigh` 的误差数值。
- 因此拒绝将该简化迭代实现替换到生产代码。它尚缺完整的退化处理与重启行为，不能凭普通随机案例通过就宣称 MATLAB 兼容。
- 本次 CI 的单测、参考导出与诊断均成功；最终仍由原始特征的两项严格对照返回失败。

## CI 配置问题与数值问题分开处理

- 已修复诊断脚本的工作目录冲突，确保读取本次交接目录中的 `run_reference.m`。
- 已修复 R2025b 不支持给 `RandStream.list` 接收返回值的问题。
- 目前失败点仍是上述两项原始特征对照，不能把脚本能够运行当成数值复现完成。
- MATLAB 与 Octave 自身的原始特征已经不一致。因此，要求同一份 Python 特征同时逐元素匹配两者不是有效的统一基准。此 PR 只记录事实，尚未更改这些门禁。

## 后续边界

### 后续进展：退化检测

[提交 41040e3 的真实 MATLAB CI](https://github.com/LevelDownRefine/WEIDA/actions/runs/45205058429)中，诊断候选已从 121/126 提升为 **125/126**，与 Windows 本机一致：

- 原候选将约 1e-30 的残差继续归一化，导致失去数值独立性并放大误差。
- 新增以机器精度和正交化前残差范数为依据的退化检测，停止无效扩展；不根据目标特征调阈值，也不改变对照容差。
- 3 个恒定类案例和分块重根案例已匹配；圆周相位重根仍失败，Linux 最大误差约 0.107607。
- 新增可手算秩一矩阵与精确不变子空间的回归测试，43 项 Python 测试通过。
- 该候选仍需要 MATLAB 起点，未实现完整重启，也未替换生产算法。生产路径的两项严格对照仍然失败。

### 后续进展：独立起始向量

[提交 12c865a 的真实 MATLAB CI](https://github.com/LevelDownRefine/WEIDA/actions/runs/45208001845)验证了新的独立生成路径：

- `matlab_start_probe.py` 按公开 dSFMT 递推、初始化和 Ziggurat 方法生成默认起点，不读取 MATLAB 特征、向量或随机样本作为生成输入，也不新增生产依赖。
- 2,048 个均匀随机样本逐位一致；2,048 个正态样本最大误差 **5.33e-15**；78 维归一化起点最大误差 **6.94e-16**。
- 使用独立起点，Windows 本机及 Linux CI 都达到 **125/126** 匹配。剩余圆周相位重根案例在 Linux 的最大误差约 **0.107733**；不可据此宣称完全等价。
- 本机从完整 BOLD 预处理开始的独立候选路径达到 **210/210**，最大误差 **6.53e-15**；该路径已加入后续 CI 诊断报告。
- 新增随机流分批调用连续性及全局 NumPy 随机状态隔离测试；本机 45 项单测及 Ruff 通过。
- dSFMT 参考代码采用 BSD 许可，已附署名与完整许可。数学来源见 [dSFMT 作者实现](https://github.com/MersenneTwister-Lab/dSFMT)及 [MathWorks 的 Ziggurat 说明](https://blogs.mathworks.com/cleve/2015/05/18/the-ziggurat-random-normal-generator/)。

重根案例的进一步隔离：MATLAB 使用 MKL 2024.1 / LAPACK 3.11.0。将捕获的 41 次矩阵向量乘法在本机 MKL 2024.2.2 重算，3,198 个元素全部逐位一致；但仅替换底层库仍不能匹配最终特征。后续应继续检查正交化、投影求解与重启的精确行为，不能简单归因于库版本。

### 尚未完成的生产兼容

最后一个圆周重根案例已进一步定位，并完成同一 MATLAB 内的一 ULP 微扰验证：起点仅改变 `1.39e-17` 即令原特征改变 `0.031535`，主子空间仍一致。详见 [重根案例根因](DEGENERATE_CASE.md)。该证据不用于跳过原始特征门禁。

- 严格保留原算法时，可使用已实现的 MATLAB 特征导出/导入路径，Python 复用其余流程；该路径的真实 MATLAB 交接测试已通过。
- 纯 Python 的原始特征兼容仍待解决，生产实现未替换；重根案例、更多独立输入及支持范围均需验证。简单改成 SciPy `eigsh`、固定全局种子或统一单向量符号都不构成修复。
- 若调整 CI 基准，应明确 MATLAB 的版本及受支持路径；Octave 应验证其实际支持的数学步骤，不能把不同的原始特征包装成 MATLAB 等价。

## 复查方法

下载上述运行的 `matlab-reference-report` artifact 后，在仓库根目录执行：

```text
uv run --locked python validation/diagnose_eigs.py <artifact目录>
```

若有 Octave，先执行以下命令，再运行同一 Python 诊断，即可把两引擎直接对照加入 JSON 报告：

```text
octave --quiet --eval "addpath('validation'); diagnose_octave_eigs('<artifact目录>')"
```

诊断报告不充当 CI 放行条件；`validation/compare.py --compare-only --strict-signal --output <artifact目录>` 继续执行原有严格门禁。

MathWorks 官方说明也明确：不同机器、版本或参数可能改变特征向量的符号，重根下还可能改变其基；默认起点使用不受全局 `rng` 控制的私有随机流。见 [eigs 文档](https://www.mathworks.com/help/matlab/ref/eigs.html)。
