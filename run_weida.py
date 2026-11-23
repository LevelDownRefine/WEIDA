"""Run WEiDA on AAL1ctx CSV files under a derivatives/fmripost directory."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from weida import (
    analyze_state_occupancy,
    calculate_state_dynamics,
    cluster_brain_states,
    compute_weighted_eigenvectors,
)


def load_scans(root, exclude_stn=()):
    scans, metadata, header = [], [], None
    root = Path(root)
    files = sorted(root.glob('sub-*/*_BOLD_AAL1ctx*.csv'))
    if (root / 'time_series').is_dir():
        files = []
        header = (root / 'rois_labels.txt').read_text(encoding='utf-8-sig').splitlines()
        if len(header) != 78 or len(set(header)) != 78:
            raise ValueError('Expected 78 unique AAL1ctx region labels')
        conditions = {'HC': 'HC', 'STNoff': 'OFF', 'STNon': 'ON'}
        with (root / 'metadata.csv').open(newline='', encoding='utf-8-sig') as handle:
            records = sorted(csv.DictReader(handle), key=lambda row: row['subject_id'])
        for row in records:
            subject = row['subject_id']
            if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn):
                continue
            path = root / 'time_series' / subject / 'tseries.csv'
            x = np.loadtxt(path, delimiter=',', ndmin=2)
            if x.shape[0] != len(header):
                raise ValueError(f'Unexpected region count: {path}')
            scans.append(x)
            metadata.append({'subject': subject, 'condition': conditions[row['condition']],
                             'file': str(path), 'frames': x.shape[1]})
    for path in files:
        subject = path.parent.name
        if subject.startswith('sub-HC'): condition = 'HC'
        elif subject.startswith('sub-STNoff'): condition = 'OFF'
        elif subject.startswith('sub-STNon'): condition = 'ON'
        else: continue
        # Exclude the invalid scan and its pair explicitly.
        if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn): continue
        with path.open(newline='', encoding='utf-8-sig') as f:
            rows = list(csv.reader(f))
        if len(rows) < 2:
            raise ValueError(f'CSV 缺少表头或时间序列数据: {path}')
        names = rows[0][1:]
        if len(names) != 78 or len(set(names)) != 78 or (header is not None and names != header):
            raise ValueError(f'Unexpected atlas columns: {path}')
        header = names
        x = np.array([[float(v) for v in row[1:]] for row in rows[1:]]).T
        if x.shape[0] != 78: raise ValueError(f'Unexpected shape: {path}')
        scans.append(x); metadata.append({'subject': subject, 'condition': condition, 'file': str(path), 'frames': x.shape[1]})
    if not scans: raise ValueError('No matching scans found')
    off = {m['subject'].removeprefix('sub-STNoff') for m in metadata if m['condition']=='OFF'}
    on = {m['subject'].removeprefix('sub-STNon') for m in metadata if m['condition']=='ON'}
    if off != on: raise ValueError(f'Unpaired STN IDs: {off ^ on}')
    ids = [m['subject'] for m in metadata]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate scans for one subject/condition')
    return scans, metadata, header


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='fmripost directory or AAL1ctx PyLEiDA input directory')
    parser.add_argument('output', help='output directory')
    parser.add_argument('--tr', type=float, default=2.1)
    parser.add_argument('--replicates', type=int, default=300)
    parser.add_argument('--seed', type=int, default=0)
    parser.add_argument('--exclude-stn', nargs='*', default=[], help='Explicit paired exclusions, e.g. 18')
    args = parser.parse_args()
    scans, meta, names = load_scans(args.input, args.exclude_stn)
    vectors, lengths = compute_weighted_eigenvectors(scans, tr=args.tr)
    clustering = cluster_brain_states(vectors, random_state=args.seed, replicates=args.replicates)
    labels, centers = clustering['labels'], clustering['centers']
    occupancy, lifetime = analyze_state_occupancy(labels, lengths, tr=args.tr)
    result = calculate_state_dynamics(labels, centers, lengths, occupancy, lifetime, [m['condition'] for m in meta])
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / 'clustering.npz', **clustering)
    np.savez_compressed(out / 'results.npz', vectors=vectors, lengths=lengths, **{k:v for k,v in result.items() if k != 'group_transitions'}, **{'group_'+k:v for k,v in result['group_transitions'].items()})
    (out / 'metadata.json').write_text(json.dumps({'scans':meta, 'regions':names, 'tr':args.tr, 'seed':args.seed, 'replicates':args.replicates, 'excluded_stn':args.exclude_stn, 'eigenvector_sign':'largest-absolute-component-positive', 'band_hz':[0.04,0.07], 'state_index_base':0}, indent=2), encoding='utf-8')

    np.savetxt(out/'Vemp.csv', result['centers'].T, delimiter=',')
    np.savetxt(out/'State_ts.csv', result['states']+1, delimiter=',', fmt='%d')
    for i, matrix in enumerate(result['transitions'],1):
        np.savetxt(out/f'Ptrans{i}.csv',matrix,delimiter=',')
    for i, condition in enumerate(['HC','OFF','ON'],1):
        selected=np.array([m['condition']==condition for m in meta])
        # Compact group rows, unlike MATLAB's all-scan NaN placeholders.
        np.savetxt(out/f'P{i}emp.csv',result['occupancy'][selected],delimiter=',')
        np.savetxt(out/f'LT{i}emp.csv',result['lifetime'][selected],delimiter=',')
        if condition in result['group_transitions']:
            np.savetxt(out/f'PTR{i}emp.csv',result['group_transitions'][condition],delimiter=',')

if __name__ == '__main__': main()
