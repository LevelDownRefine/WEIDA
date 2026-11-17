"""Four stages of the WEiDA state-space analysis."""
import numpy as np
from scipy.signal import butter, detrend, filtfilt, hilbert
from sklearn.cluster import KMeans


def compute_weighted_eigenvectors(scans, tr=2.1, band=(0.04, 0.07)):
    """Calculate framewise weighted eigenvectors from region-by-time scans."""
    b, a = butter(2, band, btype='bandpass', fs=1 / tr)
    result, lengths = [], []
    for scan in scans:
        x = np.asarray(scan, dtype=float)
        if x.ndim != 2 or not np.isfinite(x).all() or np.any(np.ptp(x, axis=1) == 0):
            raise ValueError('Scan must contain finite, nonconstant region time series')
        phase = np.angle(hilbert(filtfilt(b, a, detrend(x, axis=1), axis=1), axis=1))
        for phi in phase.T:
            matrix = np.cos(phi[:, None] - phi[None, :])
            values, vectors = np.linalg.eigh(matrix)
            idx = np.argsort(np.abs(values))[-min(6, len(values)):][::-1]
            v = vectors[:, idx] @ values[idx] / values[idx].sum()
            positive = np.mean(v > 0)
            if positive > .5 or (positive == .5 and v[v > 0].sum() > -v[v < 0].sum()):
                v = -v
            result.append(v)
        lengths.append(x.shape[1])
    return np.asarray(result), np.asarray(lengths)


def cluster_brain_states(vectors, k=4, random_state=0, replicates=300):
    """Cluster vectors using squared Euclidean distance."""
    model = KMeans(n_clusters=k, n_init=replicates, max_iter=500, random_state=random_state).fit(vectors)
    return model.labels_, model.cluster_centers_


def analyze_state_occupancy(labels, lengths, k=4, tr=2.1):
    """Per-scan occupancy and dwell time, excluding boundary runs."""
    if sum(lengths) != len(labels):
        raise ValueError('Scan lengths do not match labels')
    occupancy = np.empty((len(lengths), k)); lifetime = np.full_like(occupancy, np.nan)
    for s, seq in enumerate(np.split(labels, np.cumsum(lengths)[:-1])):
        for state in range(k):
            active = seq == state
            occupancy[s, state] = active.mean()
            diff = np.diff(active.astype(int))
            starts = np.flatnonzero(diff == 1); ends = np.flatnonzero(diff == -1)
            if len(ends) and len(starts) and ends[0] < starts[0]: ends = ends[1:]
            if len(starts) > len(ends): starts = starts[:len(ends)]
            if len(starts): lifetime[s, state] = np.mean(ends - starts) * tr
    return occupancy, lifetime


def calculate_state_dynamics(labels, centers, lengths, occupancy, lifetime, conditions):
    """Sort states by frequency and compute per-scan and group transitions."""
    k = len(centers)
    if len(conditions) != len(lengths) or sum(lengths) != len(labels):
        raise ValueError('Inconsistent scan metadata')
    order = np.argsort(-np.bincount(labels, minlength=k), kind='stable')
    remap = np.empty(k, dtype=int); remap[order] = np.arange(k)
    states = remap[labels]
    transitions = np.full((len(lengths), k, k), np.nan)
    for s, seq in enumerate(np.split(states, np.cumsum(lengths)[:-1])):
        counts = np.zeros((k, k)); np.add.at(counts, (seq[:-1], seq[1:]), 1)
        totals = counts.sum(axis=1)
        np.divide(counts, totals[:, None], out=transitions[s], where=totals[:, None] != 0)
    group = {}
    for condition in np.unique(conditions):
        matrices = transitions[np.asarray(conditions) == condition]
        valid = np.isfinite(matrices); sums = np.nansum(matrices, axis=0); counts = valid.sum(axis=0)
        group[condition] = np.divide(sums, counts, out=np.full((k,k), np.nan), where=counts != 0)
    return dict(centers=centers[order], occupancy=occupancy[:, order], lifetime=lifetime[:, order], states=states, transitions=transitions, group_transitions=group)
