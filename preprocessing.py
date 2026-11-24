"""共享 BOLD 幅值预处理，与后续分析方法无关。"""
import numpy as np
from weida_preprocessing import preprocess_scan


def preprocess_bold(scan, tr=2.1, band=(0.04, 0.07)):
    """返回脑区 × 时间的幅值矩阵，不修改输入。

    按原版 WEiDA 顺序执行：去线性趋势、去均值、二阶 Butterworth 带通，
    再用 filtfilt 双向滤波；默认滤波器的 padlen 为 12。
    填充仅用于临时延拓，不删除原始时间点；此处不提取相位、不裁剪首尾，
    也不缓存到磁盘。band=None 显式保留输入幅值，不去趋势、不滤波。
    """
    x = np.array(scan, dtype=float, copy=True)
    if x.ndim != 2 or min(x.shape) == 0 or not np.isfinite(x).all():
        raise ValueError('Expected nonempty finite region x time array')
    if not np.isfinite(tr) or tr <= 0 or np.any(np.ptp(x, axis=1) == 0):
        raise ValueError('Invalid TR or constant region')
    if band is None:
        return x
    if not 0 < band[0] < band[1] < .5/tr:
        raise ValueError('Invalid bandpass')
    return preprocess_scan(x, remove_trend=True, demean=True, band=band, tr=tr)
