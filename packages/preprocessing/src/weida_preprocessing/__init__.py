"""BOLD 幅值预处理；不提取相位，不裁剪首尾时间点。"""
import numpy as np
from scipy.signal import butter, detrend, filtfilt


def preprocess_scan(scan, *, remove_trend=False, demean=False, band=None, tr=None):
    """处理一个脑区 × 时间矩阵；各项幅值变换均需显式启用。"""
    data = np.array(scan, dtype=float, copy=True)
    if data.ndim != 2 or not data.size or not np.isfinite(data).all():
        raise ValueError('Expected a finite region-by-time matrix')
    if remove_trend:
        data = detrend(data, axis=1)
    if demean:
        data -= data.mean(axis=1, keepdims=True)
    if band is not None:
        if tr is None or not np.isfinite(tr) or tr <= 0:
            raise ValueError('Bandpass requires a positive TR in seconds')
        low, high = band
        if not 0 < low < high < 0.5 / tr:
            raise ValueError('Band must lie between zero and the Nyquist frequency')
        b, a = butter(2, [low, high], btype='bandpass', fs=1 / tr)
        padlen = 3 * (max(len(a), len(b)) - 1)
        if data.shape[1] <= padlen:
            raise ValueError("Insufficient time points")
        data = filtfilt(b, a, data, axis=1, padlen=padlen)
    return data
