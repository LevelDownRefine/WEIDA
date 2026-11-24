# MATLAB / Octave 委托测试说明

此目录只有合成数据，不包含受试者数据。不要修改 input.mat 或测试脚本。

1. 解压整个目录；在 MATLAB 中将当前目录切换到这里。
2. 需要 Signal Processing Toolbox 和 Statistics and Machine Learning Toolbox；建议 MATLAB R2020b 或更新版本。
3. 在命令窗口运行 `run_reference`。不需要 Python，不需要原始影像或完整研究数据。
4. 把整个目录压缩回传，至少保留 input.mat、manifest.json、reference.m、run_reference.m、reference.mat 和 runtime.log。
5. 如果报错，不要自行改算法；把 runtime.log 和完整报错回传。

命令行也可使用 `matlab -batch "run_reference"`。当前交接与验收使用真实 MATLAB；Octave 不再作为验证基准。

本任务是独立四阶段参考实现的对照，不是运行整份原始论文脚本。保留 MATLAB 求解器返回的特征向量符号，仅按原脚本调整加权后总向量符号；原始特征结果与 Python 的对照为必需检查。

接收方在 Python 仓库运行：

```
uv run --locked python validation/compare.py --compare-only --output /path/to/returned-folder
```

该命令先验证输入和脚本的 SHA256，再生成 report.json。不得把 MATLAB 对照通过宣称为已复现论文全部结果。
