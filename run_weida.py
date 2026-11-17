"""Run WEiDA on AAL1ctx CSV files under a derivatives/fmripost directory."""
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from weida import compute_weighted_eigenvectors, cluster_brain_states, analyze_state_occupancy, calculate_state_dynamics


def load_scans(root):
    scans, metadata, header = [], [], None
    files = sorted(Path(root).glob('sub-*/*_BOLD_AAL1ctx*.csv'))
    for path in files:
        subject = path.parent.name
        if subject.startswith('sub-HC'): condition = 'HC'
        elif subject.startswith('sub-STNoff'): condition = 'OFF'
        elif subject.startswith('sub-STNon'): condition = 'ON'
        else: continue
        # Exclude the invalid scan and its pair explicitly.
        if subject in ('sub-STNon18', 'sub-STNoff18'): continue
        with path.open(newline='', encoding='utf-8-sig') as f:
            rows = list(csv.reader(f))
        names = rows[0][1:]
        if len(names) != 78 or (header is not None and names != header):
            raise ValueError(f'Unexpected atlas columns: {path}')
        header = names
        x = np.array([[float(v) for v in row[1:]] for row in rows[1:]]).T
        if x.shape[0] != 78: raise ValueError(f'Unexpected shape: {path}')
        scans.append(x); metadata.append({'subject': subject, 'condition': condition, 'file': str(path), 'frames': x.shape[1]})
    if not scans: raise ValueError('No matching scans found')
    return scans, metadata, header


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='fmripost directory')
    parser.add_argument('output', help='output directory')
    parser.add_argument('--tr', type=float, default=2.1)
    parser.add_argument('--replicates', type=int, default=300)
    parser.add_argument('--seed', type=int, default=0)
    args = parser.parse_args()
    scans, meta, names = load_scans(args.input)
    vectors, lengths = compute_weighted_eigenvectors(scans, tr=args.tr)
    labels, centers = cluster_brain_states(vectors, random_state=args.seed, replicates=args.replicates)
    occupancy, lifetime = analyze_state_occupancy(labels, lengths, tr=args.tr)
    result = calculate_state_dynamics(labels, centers, lengths, occupancy, lifetime, [m['condition'] for m in meta])
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / 'results.npz', vectors=vectors, lengths=lengths, **{k:v for k,v in result.items() if k != 'group_transitions'}, **{'group_'+k:v for k,v in result['group_transitions'].items()})
    (out / 'metadata.json').write_text(json.dumps({'scans':meta, 'regions':names, 'tr':args.tr, 'seed':args.seed, 'replicates':args.replicates}, indent=2), encoding='utf-8')

if __name__ == '__main__': main()
