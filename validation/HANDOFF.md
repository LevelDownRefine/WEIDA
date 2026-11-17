# MATLAB / Octave 委托测试说明

此目录只有合成数据，不包含受试者数据。不要修改 input.mat 或测试脚本。

1. 解压整个目录；在 MATLAB 中将当前目录切换到这里。
2. 需要 Signal Processing Toolbox 和 Statistics and Machine Learning Toolbox；建议 MATLAB R2020b 或更新版本。
3. 在命令窗口运行 `run_reference`。不需要 Python，不需要原始影像或完整研究数据。
4. 把整个目录压缩回传，至少保留 input.mat、manifest.json、reference.m、run_reference.m、reference.mat 和 runtime.log。
5. 如果报错，不要自行改算法；把 runtime.log 和完整报错回传。

命令行也可使用 `matlab -batch "run_reference"`。Octave 使用 `octave --no-gui --quiet --eval "run_reference"`，需 signal/statistics 包。

本任务是独立四阶段参考实现的对照，不是运行整份原始论文脚本。MATLAB 内置函数的原始特征向量符号结果与统一符号后的结果会分别导出；不要求前者与 Python 完全一致。

接收方在 Python 仓库运行：

```
uv run --locked python validation/compare.py --compare-only --output /path/to/returned-folder
```

该命令先验证输入和脚本的 SHA256，再生成 report.json。不得把 MATLAB 对照通过宣称为已复现论文全部结果。
