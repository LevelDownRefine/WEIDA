"""Amplitude preprocessing shared by algorithms that opt into this profile."""
import numpy as np
from scipy.signal import butter, detrend, filtfilt


def preprocess_bold(scan, tr=2.1, band=(0.04, 0.07)):
    """Return region × time amplitudes; never convert amplitudes to phase.

    Linear detrend, demean, second-order Butterworth bandpass and zero-phase
    filtering with MATLAB-style pad length. Input is not modified.
    """
    band = np.asarray(band, dtype=float)
    if (not np.isfinite(tr) or tr <= 0 or band.shape != (2,)
            or not np.isfinite(band).all() or not 0 < band[0] < band[1] < .5 / tr):
        raise ValueError('Invalid TR or bandpass')
    b, a = butter(2, band, btype='bandpass', fs=1 / tr)
    padlen = 3 * (max(len(a), len(b)) - 1)
    x = np.asarray(scan, dtype=float)
    if x.ndim != 2 or min(x.shape) == 0 or not np.isfinite(x).all():
        raise ValueError('Expected nonempty finite region x time array')
    if np.any(np.ptp(x, axis=1) == 0) or x.shape[1] <= padlen:
        raise ValueError('Constant region or insufficient time points')
    ts = detrend(x, axis=1)
    ts -= ts.mean(axis=1, keepdims=True)
    return filtfilt(b, a, ts, axis=1, padlen=padlen)
