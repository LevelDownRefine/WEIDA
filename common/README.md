# 共用前后处理

- `preprocessing/`：BOLD 幅值预处理，包名 `weida-preprocessing`。
- `postprocessing/`：指标配对、组间统计、校正、通用作图及保存，包名 `weida-postprocessing`。

两个目录分别可安装；调用方固定 Git 提交及对应子目录，只安装需要的包。Python 导入名保持 `weida_preprocessing` 和 `weida_postprocessing`。测试统一放在根目录 `tests/`，算法专属实现不放在这里。
