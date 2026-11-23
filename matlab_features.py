"""Prepare shared amplitudes for MATLAB and validate returned WEiDA features."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.io import loadmat, savemat

from preprocessing import preprocess_bold


def input_signature(scans, metadata, regions, tr):
    """Bind features to scan/ROI order, raw amplitudes and preprocessing settings."""
    info = {'subjects': [m['subject'] for m in metadata],
            'conditions': [m['condition'] for m in metadata], 'regions': list(regions),
            'shapes': [list(x.shape) for x in scans], 'tr': float(tr),
            'preprocessing': 'linear-detrend-demean-butter2-filtfilt-pad12', 'band': [.04, .07]}
    digest = hashlib.sha256(json.dumps(info, sort_keys=True).encode('utf-8'))
    for scan in scans:
        digest.update(np.asarray(scan, dtype='<f8').tobytes(order='C'))
    return digest.hexdigest()


def prepare(scans, metadata, regions, tr, output):
    cells = np.empty((len(scans), 1), dtype=object)
    for i, scan in enumerate(scans):
        cells[i, 0] = preprocess_bold(scan, tr=tr)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    savemat(output / 'input.mat', {
        'scans': cells, 'tr': float(tr),
        'subjects': np.asarray([m['subject'] for m in metadata], dtype=object),
        'regions': np.asarray(regions, dtype=object),
        'input_signature': input_signature(scans, metadata, regions, tr),
    })
    (output / 'export_weida_features.m').write_bytes(
        Path(__file__).with_name('export_weida_features.m').read_bytes())


def load_features(path, scans, metadata, regions, tr):
    data = loadmat(path, simplify_cells=True)
    if str(data['input_signature']) != input_signature(scans, metadata, regions, tr):
        raise ValueError('MATLAB features do not match the input scans, ROI order or TR')
    subjects = np.atleast_1d(data['subjects']).tolist()
    labels = np.atleast_1d(data['regions']).tolist()
    if subjects != [m['subject'] for m in metadata] or labels != list(regions):
        raise ValueError('MATLAB feature subject/ROI order mismatch')
    lengths = np.atleast_1d(data['lengths'])
    expected = np.asarray([x.shape[1] for x in scans])
    if not np.array_equal(lengths, expected) or float(data['tr']) != tr:
        raise ValueError('MATLAB feature frame count or TR mismatch')
    # Do not squeeze the feature axes, even for a one-frame input.
    vectors = loadmat(path, variable_names=['vectors'])['vectors']
    if (vectors.shape != (int(expected.sum()), len(regions)) or np.iscomplexobj(vectors)
            or not np.isfinite(vectors).all()):
        raise ValueError('Invalid MATLAB feature matrix')
    engine = str(data['runtime_engine'])
    if engine not in ('MATLAB', 'Octave'):
        raise ValueError('Unknown feature runtime')
    return vectors, expected, {'engine': engine, 'version': str(data['runtime_version']),
                               'file': str(Path(path).resolve()), 'input_signature': str(data['input_signature'])}


def main():
    from run_weida import load_scans

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input')
    parser.add_argument('output', help='New directory for the MATLAB handoff')
    parser.add_argument('--tr', type=float, default=2.1)
    parser.add_argument('--exclude-stn', nargs='*', default=[])
    args = parser.parse_args()
    scans, metadata, regions = load_scans(args.input, args.exclude_stn)
    prepare(scans, metadata, regions, args.tr, args.output)
    print(f'Handoff ready: {Path(args.output).resolve()}')


if __name__ == '__main__':
    main()
