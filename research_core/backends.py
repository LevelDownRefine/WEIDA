"""Adapters: same input/output contract, distinct WEiDA and official PhiID algorithms."""
import numpy as np
from .contracts import AlgorithmResult

ATOM_NAMES = tuple(a + 't' + b for a in 'rxys' for b in 'rxys')
SWAP_XY = np.array([ATOM_NAMES.index(a.translate(str.maketrans('xy', 'yx')))
                    for a in ATOM_NAMES])
PHIID_REVISION = '6c5f2e9d33c985efbdf875d45cb5a2a6a5cdbf44'


def weida_backend(scans, metadata, regions, *, tr, k=4, seed=0, replicates=300):
    from weida import (weighted_eigenvectors_from_amplitudes, cluster_brain_states,
                       analyze_state_occupancy, calculate_state_dynamics)
    from downstream import state_metastability

    if any(isinstance(x, bool) or not isinstance(x, (int, np.integer)) or x < 1
           for x in (k, replicates)) or k < 2 or k >= sum(x.shape[1] for x in scans):
        raise ValueError('Require 2 <= k < total frames and positive integer replicates')
    vectors, lengths = weighted_eigenvectors_from_amplitudes(scans)
    clustering = cluster_brain_states(vectors, k=k, random_state=seed, replicates=replicates)
    occupancy, lifetime = analyze_state_occupancy(clustering['labels'], lengths, k=k, tr=tr)
    result = calculate_state_dynamics(clustering['labels'], clustering['centers'], lengths,
                                     occupancy, lifetime, [m['condition'] for m in metadata])
    transitions = np.where(result['transition_valid'][:, :, None], result['transitions'], np.nan)
    features = dict(occupancy=result['occupancy'], lifetime=result['lifetime'],
                    transitions=transitions.reshape(len(scans), -1),
                    metastability=state_metastability(vectors, result['states'], result['centers'], lengths))
    states = [f'state_{i+1}' for i in range(k)]
    artifacts = {key: value for key, value in result.items() if key != 'group_transitions'}
    artifacts.update(vectors=vectors, lengths=lengths,
                     **{'group_' + key: value for key, value in result['group_transitions'].items()})
    return AlgorithmResult('weida', features,
        dict(occupancy=states, lifetime=states, metastability=states,
             transitions=[f'{i+1}_to_{j+1}' for i in range(k) for j in range(k)]),
        dict(occupancy='fraction', lifetime='s', transitions='probability', metastability='WEiDA variance'),
        dict(occupancy='Fraction of scan frames in each state',
             lifetime='Mean complete-run dwell time; boundary runs excluded',
             transitions='Per-scan transitions; unobserved source rows excluded as NaN',
             metastability='Squared distance to state center / (regions * (state frames - 1))'),
        artifacts, dict(k=k, seed=seed, replicates=replicates,
                        eigenvector_sign='largest-absolute-component-positive'))


def phiid_backend(scans, metadata, regions, *, tr, lag=1):
    # Official estimator only; no local entropy/redundancy decomposition implementation.
    from phyid.calculate import calc_PhiID

    if isinstance(lag, bool) or not isinstance(lag, (int, np.integer)) or lag < 1:
        raise ValueError('lag must be a positive integer in frames')
    left, right = np.triu_indices(len(regions), 1)
    upper = np.empty((len(scans), len(left), 16))
    for s, scan in enumerate(scans):
        if scan.shape[1] - lag <= 4:
            raise ValueError(f"{metadata[s]['subject']}: at least five lagged pairs required")
        for edge, (i, j) in enumerate(zip(left, right)):
            data = scan[[i, j]]
            paired = np.vstack((data[:, :-lag], data[:, lag:]))
            if (np.any(np.std(paired, axis=1) == 0)
                    or np.linalg.eigvalsh(np.corrcoef(paired)).min() <= 1e-12):
                raise ValueError(f"{metadata[s]['subject']}, {regions[i]}/{regions[j]}: singular lagged covariance")
            atoms, _ = calc_PhiID(data[0], data[1], int(lag), kind='gaussian', redundancy='MMI')
            local = np.column_stack([atoms[name] for name in ATOM_NAMES]) / np.log(2)
            if not np.isfinite(local).all():
                raise ValueError(f"{metadata[s]['subject']}: upstream returned nonfinite atoms")
            upper[s, edge] = local.mean(axis=0)
    matrices = np.full((len(scans), len(regions), len(regions), 16), np.nan)
    matrices[:, left, right, :] = upper
    matrices[:, right, left, :] = upper[..., SWAP_XY]
    # Both endpoint assignments contribute equally; never average only upper triangles.
    mean = np.concatenate((upper, upper[..., SWAP_XY]), axis=1).mean(axis=1)
    return AlgorithmResult('phiid', {'phiid_mean': mean}, {'phiid_mean': list(ATOM_NAMES)},
        {'phiid_mean': 'bit'},
        {'phiid_mean': 'Gaussian MMI; mean over all nonself edges with both endpoint assignments'},
        dict(atoms=upper, matrices=matrices, roi_x=left, roi_y=right, atom_names=np.array(ATOM_NAMES)),
        dict(lag=int(lag), lag_seconds=float(lag*tr), kind='gaussian', redundancy='MMI',
             upstream='Imperial-MIND-lab/integrated-info-decomp', upstream_revision=PHIID_REVISION))


BACKENDS = {'weida': weida_backend, 'phiid': phiid_backend}
