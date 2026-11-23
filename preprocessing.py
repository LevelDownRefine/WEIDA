"""Shared BOLD amplitude preprocessing, independent of the analysis method."""
import numpy as np
from scipy.signal import butter, detrend, filtfilt


def preprocess_bold(scan, tr=2.1, band=(0.04, 0.07)):
    """Return region x time amplitudes without modifying the input.

    Match the original WEiDA order: linear detrend, demean, second-order
    Butterworth bandpass, then filtfilt with padlen=12 for the default filter.
    Padding is temporary extension, not removal of original frames. No phase
    extraction, edge trimming or disk caching is performed here.
    band=None explicitly preserves loaded amplitudes (no detrend or filtering).
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
    b, a = butter(2, band, btype='bandpass', fs=1/tr)
    padlen = 3*(max(len(a), len(b))-1)
    if x.shape[1] <= padlen:
        raise ValueError('Insufficient time points')
    ts = detrend(x, axis=1)
    ts -= ts.mean(axis=1, keepdims=True)
    return filtfilt(b, a, ts, axis=1, padlen=padlen)
