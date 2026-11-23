"""同一读取/预处理流程，仅切换计算方法；适用的统计和图形共用。"""
import argparse
import json
from pathlib import Path

import numpy as np

from downstream import independent_comparison, paired_comparison, state_metastability
from run_weida import load_scans
from weida import (
    analyze_state_occupancy,
    calculate_state_dynamics,
    cluster_brain_states,
    compute_weighted_eigenvectors,
    preprocess_bold,
)


def weida_method(scans, conditions, tr, *, k=4, seed=0, replicates=300):
    """Accept processed amplitudes; keep phase and state modelling inside WEiDA."""
    vectors, lengths = compute_weighted_eigenvectors(scans, tr=tr, preprocessed=True)
    clustering = cluster_brain_states(vectors, k=k, random_state=seed, replicates=replicates)
    p, lt = analyze_state_occupancy(clustering['labels'], lengths, k=k, tr=tr)
    native = calculate_state_dynamics(clustering['labels'], clustering['centers'], lengths, p, lt, conditions)
    groups = native.pop('group_transitions')
    native.update(vectors=vectors, lengths=lengths, **{'group_'+g: v for g, v in groups.items()})
    trans = np.where(native['transition_valid'][:, :, None], native['transitions'], np.nan)
    names = [f'state_{i+1}' for i in range(k)]
    # A family is (scan x feature values, column labels, physical unit).
    features = {
        'occupancy': (native['occupancy'], names, 'fraction'),
        'lifetime': (native['lifetime'], names, 's'),
        'transitions': (trans.reshape(len(scans), -1), [f'{i+1}_to_{j+1}' for i in range(k) for j in range(k)], 'probability'),
        'metastability': (state_metastability(vectors, native['states'], native['centers'], lengths), names, 'WEiDA variance'),
    }
    return native, features, {'k': k, 'seed': seed, 'replicates': replicates}


def phiid_method(scans, conditions, tr, *, lag=1):
    """Call official Gaussian MMI; this adapter contains no PhiID decomposition formula."""
    from phyid.calculate import calc_PhiID

    if isinstance(lag, bool) or not isinstance(lag, (int, np.integer)) or lag < 1:
        raise ValueError('lag must be a positive integer')
    names = tuple(a+'t'+b for a in 'rxys' for b in 'rxys')
    swap = [names.index(a.translate(str.maketrans('xy', 'yx'))) for a in names]
    n = scans[0].shape[0]
    left, right = np.triu_indices(n, 1)
    atoms = np.empty((len(scans), len(left), 16))
    for s, scan in enumerate(scans):
        if scan.shape[1]-lag <= 4:
            raise ValueError(f'Scan {s}: at least five lagged pairs required')
        for edge, (i, j) in enumerate(zip(left, right)):
            data = scan[[i, j]]
            paired = np.vstack((data[:, :-lag], data[:, lag:]))
            if np.any(np.std(paired, axis=1) == 0) or np.linalg.eigvalsh(np.corrcoef(paired)).min() <= 1e-12:
                raise ValueError(f'Scan {s}, regions {i}/{j}: singular lagged covariance')
            local, _ = calc_PhiID(data[0], data[1], int(lag), kind='gaussian', redundancy='MMI')
            bits = np.column_stack([local[name] for name in names]) / np.log(2)
            if not np.isfinite(bits).all():
                raise ValueError('Official PhiID returned nonfinite values')
            atoms[s, edge] = bits.mean(axis=0)
    matrices = np.full((len(scans), n, n, 16), np.nan)
    matrices[:, left, right] = atoms
    matrices[:, right, left] = atoms[..., swap]
    means = np.concatenate((atoms, atoms[..., swap]), axis=1).mean(axis=1)
    return ({'atoms': atoms, 'matrices': matrices, 'roi_x': left, 'roi_y': right, 'atom_names': np.array(names)},
            {'phiid_mean': (means, names, 'bit')},
            {'lag': int(lag), 'lag_seconds': float(lag*tr), 'kind': 'gaussian', 'redundancy': 'MMI',
             'upstream_revision': '6c5f2e9d33c985efbdf875d45cb5a2a6a5cdbf44'})


def cohort_indices(metadata):
    groups = {g: {} for g in ('HC', 'OFF', 'ON')}
    prefixes = {'HC': 'sub-HC', 'OFF': 'sub-STNoff', 'ON': 'sub-STNon'}
    ids = []
    for i, row in enumerate(metadata):
        group = row['condition']
        if group not in groups or not row['subject'].startswith(prefixes[group]):
            raise ValueError('Expected HC/STN subject names and HC/OFF/ON conditions')
        subject = row['subject'][len(prefixes[group]):]
        if not subject or subject in groups[group]:
            raise ValueError('Empty or duplicate participant-condition')
        groups[group][subject] = i
        ids.append(subject)
    if groups['OFF'].keys() != groups['ON'].keys() or any(len(g) < 2 for g in groups.values()):
        raise ValueError('Require complete ON/OFF pairs and at least two participants per group')
    return ids, {g: [v[s] for s in sorted(v)] for g, v in groups.items()}


def analyze(scans, metadata, *, method, preprocessing='weida_amplitude', tr=2.1, **options):
    """Only the middle calculation varies; preprocessing and comparisons are shared."""
    methods = {'weida': weida_method, 'phiid': phiid_method}
    if method not in methods or preprocessing not in ('loaded', 'weida_amplitude'):
        raise ValueError('Unknown method or preprocessing profile')
    if len(scans) != len(metadata) or not scans:
        raise ValueError('Scan metadata mismatch')
    ids, _ = cohort_indices(metadata)
    conditions = [m['condition'] for m in metadata]
    processed = [preprocess_bold(x, tr=tr, band=(.04, .07) if preprocessing == 'weida_amplitude' else None)
                 for x in scans]
    if processed[0].shape[0] < 2 or any(x.shape[0] != processed[0].shape[0] for x in processed):
        raise ValueError('Inconsistent region counts')
    native, features, parameters = methods[method](processed, conditions, tr, **options)
    comparisons = {}
    for family, (values, _, _) in features.items():
        comparisons[family] = {
            'ON_vs_OFF': paired_comparison(values, ids, conditions),
            **{g+'_vs_HC': independent_comparison(values, conditions, target=g) for g in ('OFF', 'ON')},
        }
    return native, features, comparisons, parameters


def plot_features(features, metadata, output):
    """Shared descriptive plots only; state and connection maps remain method-specific."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    _, groups = cohort_indices(metadata)
    for family, (values, names, unit) in features.items():
        for start in range(0, len(names), 16):
            columns = range(start, min(start+16, len(names)))
            for paired in (False, True):
                fig, axes = plt.subplots((len(columns)+3)//4, min(4, len(columns)), squeeze=False, figsize=(12, 3*((len(columns)+3)//4)))
                for j, ax in zip(columns, axes.flat):
                    if paired:
                        for a, b in zip(values[groups['OFF'], j], values[groups['ON'], j]):
                            if np.isfinite([a, b]).all():
                                ax.plot([0, 1], [a, b], alpha=.5, lw=.8)
                        ax.set_xticks([0, 1], ['OFF', 'ON'])
                    else:
                        for pos, group in enumerate(('HC', 'OFF', 'ON')):
                            y = values[groups[group], j]
                            ax.scatter(np.full(len(y), pos), y, s=10)
                        ax.set_xticks([0, 1, 2], ['HC', 'OFF', 'ON'])
                    ax.set_title(names[j]); ax.set_ylabel(unit)
                for ax in list(axes.flat)[len(columns):]:
                    ax.set_visible(False)
                fig.suptitle(f'{family} | descriptive'); fig.tight_layout()
                for ext in ('png', 'pdf'):
                    fig.savefig(output / f'{family}_{start}_{"paired" if paired else "groups"}.{ext}', dpi=160)
                plt.close(fig)


def json_safe(value):
    if isinstance(value, dict):
        return {k: json_safe(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    return None if isinstance(value, float) and not np.isfinite(value) else value


def run_study(input_path, output, *, method, preprocessing='weida_amplitude', tr=2.1, exclude_stn=(), plots=True, **options):
    output = Path(output)
    if output.exists():
        raise FileExistsError(output)
    scans, metadata, regions = load_scans(input_path, exclude_stn)
    native, features, comparisons, parameters = analyze(scans, metadata, method=method, preprocessing=preprocessing, tr=tr, **options)
    output.mkdir(parents=True)
    np.savez_compressed(output/'results.npz', **native)
    np.savez_compressed(output/'features.npz', **{k: v[0] for k, v in features.items()})
    info = {'method': method, 'scans': metadata, 'regions': regions, 'tr': tr, 'preprocessing': preprocessing,
            'band_hz': [.04, .07] if preprocessing == 'weida_amplitude' else None,
            'excluded_stn': list(exclude_stn), 'parameters': parameters,
            'features': {k: {'names': v[1], 'unit': v[2]} for k, v in features.items()},
            'inference': 'Two-sided paired/equal-variance independent t; BH over finite p values per contrast and family',
            'phiid_summary': 'Mean of all nonself edges with both endpoint assignments; not edgewise inference' if method == 'phiid' else None}
    for filename, data in [('metadata.json', info), ('statistics.json', comparisons)]:
        (output/filename).write_text(json.dumps(json_safe(data), ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    if plots:
        plot_features(features, metadata, output)
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input'); parser.add_argument('output')
    parser.add_argument('--method', required=True, choices=['weida', 'phiid'])
    parser.add_argument('--preprocessing', choices=['loaded', 'weida_amplitude'], default='weida_amplitude')
    parser.add_argument('--tr', type=float, default=2.1)
    parser.add_argument('--exclude-stn', nargs='*', default=[])
    parser.add_argument('--no-plots', action='store_true')
    for option in ('k', 'seed', 'replicates', 'lag'):
        parser.add_argument('--'+option, type=int)
    args = parser.parse_args()
    options = {name: getattr(args, name) for name in ('k', 'seed', 'replicates', 'lag') if getattr(args, name) is not None}
    allowed = {'weida': {'k', 'seed', 'replicates'}, 'phiid': {'lag'}}
    if set(options)-allowed[args.method]:
        parser.error('存在不适用于当前方法的参数')
    run_study(args.input, args.output, method=args.method, preprocessing=args.preprocessing, tr=args.tr,
              exclude_stn=args.exclude_stn, plots=not args.no_plots, **options)


if __name__ == '__main__':
    main()
