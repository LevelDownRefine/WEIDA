"""Four numbered stages from WEIDA-DBS StateSpace_modelling.m.

This implements the equations, not bitwise MATLAB eigs/kmeans equivalence.
"""
import numpy as np
from scipy.signal import hilbert
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples

from preprocessing import preprocess_bold


def compute_weighted_eigenvectors(scans, tr=2.1, band=(0.04, 0.07), *, preprocessed=False):
    """Stage 1. Return frame x region features and per-scan lengths.

    Use six largest-magnitude eigenpairs. Canonicalize EACH eigenvector by
    making its largest absolute component positive before weighting. This is
    an explicit reproducibility convention absent from the MATLAB original.
    It changes the original solver-dependent convention; validate scientifically.
    MATLAB-style filtfilt padding is 3*(max(len(a),len(b))-1).
    preprocessed=True accepts amplitudes from the shared preprocessing step.
    """
    if not np.isfinite(tr) or tr <= 0 or (not preprocessed and not 0 < band[0] < band[1] < .5/tr):
        raise ValueError('Invalid TR or bandpass')
    chunks, lengths, regions = [], [], None
    for scan in scans:
        x = np.asarray(scan, dtype=float)
        if x.ndim != 2 or min(x.shape) == 0 or not np.isfinite(x).all():
            raise ValueError('Expected nonempty finite region x time array')
        if regions is not None and regions != x.shape[0]:
            raise ValueError('Inconsistent region counts')
        regions = x.shape[0]
        ts = x if preprocessed else preprocess_bold(x, tr=tr, band=band)
        phase = np.angle(hilbert(ts, axis=1))
        chunk = []
        for phi in phase.T:
            values, vectors = np.linalg.eigh(np.cos(phi[:,None]-phi[None,:]))
            idx = np.argsort(np.abs(values))[-min(6,regions):][::-1]
            values, vectors = values[idx], vectors[:,idx]
            anchors = vectors[np.argmax(np.abs(vectors),axis=0),np.arange(len(idx))]
            vectors *= np.where(anchors < 0,-1.,1.)
            v = vectors @ values / values.sum()
            positive = np.mean(v > 0)
            if positive > .5 or (positive == .5 and v[v>0].sum() > -v[v<0].sum()): v = -v
            chunk.append(v)
        chunks.append(np.asarray(chunk)); lengths.append(x.shape[1])
    if not chunks: raise ValueError('No scans')
    return np.concatenate(chunks), np.asarray(lengths)


def cluster_brain_states(vectors, k=4, random_state=0, replicates=300):
    """Stage 2. MATLAB-like IDX/C/SUMD/D and exact squared-distance silhouette.

    Labels are zero-based internally. K-means initialization is not MATLAB-identical.
    """
    model = KMeans(n_clusters=k,n_init=replicates,max_iter=500,random_state=random_state,algorithm='lloyd').fit(vectors)
    distances = model.transform(vectors)**2
    labels = model.labels_
    silhouette = silhouette_samples(vectors,labels,metric='sqeuclidean')
    return {'labels': labels,'centers': model.cluster_centers_,'distances': distances,
                'sumd': np.bincount(labels,weights=distances[np.arange(len(labels)),labels],minlength=k),
                'silhouette': silhouette,'silhouette_mean': float(silhouette.mean())}


def _sequences(labels,lengths,k):
    labels, lengths = np.asarray(labels), np.asarray(lengths)
    if labels.ndim!=1 or lengths.ndim!=1 or not np.issubdtype(lengths.dtype,np.integer) or np.any(lengths<=0) or lengths.sum()!=len(labels):
        raise ValueError('Invalid scan lengths')
    if not np.issubdtype(labels.dtype,np.integer) or np.any((labels<0)|(labels>=k)):
        raise ValueError('Invalid state labels')
    return np.split(labels,np.cumsum(lengths)[:-1])


def analyze_state_occupancy(labels,lengths,k=4,tr=2.1):
    """Stage 3. Occupancy and complete-run dwell time in seconds."""
    seqs=_sequences(labels,lengths,k)
    occupancy=np.empty((len(seqs),k)); lifetime=np.full_like(occupancy,np.nan)
    for s,seq in enumerate(seqs):
        for c in range(k):
            active=seq==c; occupancy[s,c]=active.mean()
            changes=np.diff(active.astype(int))
            a=np.flatnonzero(changes==1); b=np.flatnonzero(changes==-1)
            if len(b)>len(a): b=b[1:]
            elif len(a)>len(b): a=a[:-1]
            elif len(a) and len(b) and a[0]>b[0]: b=b[1:]; a=a[:-1]
            if len(a) and len(b): lifetime[s,c]=np.mean(b-a)*tr
    return occupancy,lifetime


def calculate_state_dynamics(labels,centers,lengths,occupancy,lifetime,conditions):
    """Stage 4. Frequency ordering, scan transitions and eligible-scan group means.

    Missing source states have zero rows per scan (as in MATLAB). Group means
    exclude these rows; if no scan contributes, the group row is NaN.
    """
    k=len(centers); seqs=_sequences(labels,lengths,k); conditions=np.asarray(conditions)
    if conditions.shape!=(len(seqs),) or np.shape(occupancy)!=(len(seqs),k) or np.shape(lifetime)!=(len(seqs),k):
        raise ValueError('Inconsistent scan metadata')
    order=np.argsort(-np.bincount(labels,minlength=k),kind='stable')
    remap=np.empty(k,dtype=int); remap[order]=np.arange(k)
    transitions=np.zeros((len(seqs),k,k)); valid=np.zeros((len(seqs),k),dtype=bool)
    for s,seq in enumerate(seqs):
        seq=remap[seq]; counts=np.zeros((k,k)); np.add.at(counts,(seq[:-1],seq[1:]),1)
        totals=counts.sum(axis=1); valid[s]=totals>0
        np.divide(counts,totals[:,None],out=transitions[s],where=totals[:,None]>0)
    group={}
    for condition in np.unique(conditions):
        selected=conditions==condition; count=valid[selected].sum(axis=0)
        group[str(condition)]=np.divide(transitions[selected].sum(axis=0),count[:,None],out=np.full((k,k),np.nan),where=count[:,None]>0)
    return {'centers': np.asarray(centers)[order],'occupancy': occupancy[:,order],'lifetime': lifetime[:,order],'states': remap[labels],'transitions': transitions,'group_transitions': group,'order': order,'transition_valid': valid}
