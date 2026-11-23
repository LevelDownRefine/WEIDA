"""MATLAB feature handoff validation; shared fixtures also run in real MATLAB CI."""
import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.io import loadmat, savemat

from matlab_features import input_signature, load_features, prepare
from preprocessing import preprocess_bold


def make_fixture():
    rng = np.random.default_rng(753)
    scans = [rng.normal(size=(78, n)) for n in (40, 50, 45)]
    metadata = [{'subject': subject, 'condition': condition, 'frames': x.shape[1]}
                for subject, condition, x in zip(
                    ['sub-HC01', 'sub-STNoff01', 'sub-STNon01'], ['HC', 'OFF', 'ON'], scans)]
    return scans, metadata, [f'ROI{i}' for i in range(78)]


def check_export(folder, engine):
    """Called after actual MATLAB/Octave export; run Python stages 2–4 unchanged."""
    from run_weida import main

    root = Path(folder)
    scans, metadata, regions = make_fixture()
    path = root / 'features.mat'
    vectors, lengths, source = load_features(path, scans, metadata, regions, 2.1)
    assert source['engine'] == engine and source['version']
    input_root = root / 'input'
    for scan, record in zip(scans, metadata):
        directory = input_root / record['subject']
        directory.mkdir(parents=True)
        with (directory / 'fixture_BOLD_AAL1ctx.csv').open('w', newline='') as handle:
            writer = csv.writer(handle)
            writer.writerow([''] + regions)
            for i, row in enumerate(scan.T):
                writer.writerow([i] + row.tolist())
    with (patch('run_weida.compute_weighted_eigenvectors', side_effect=AssertionError('Recomputed imported features')),
          patch('sys.argv', ['run_weida.py', str(input_root), str(root / 'result'),
                             '--matlab-features', str(path), '--replicates', '2'])):
        main()
    with np.load(root / 'result' / 'results.npz') as result:
        np.testing.assert_array_equal(result['vectors'], vectors)
        np.testing.assert_array_equal(result['lengths'], lengths)
        assert result['occupancy'].shape == (3, 4)
        assert result['transitions'].shape == (3, 4, 4)
    info = json.loads((root / 'result' / 'metadata.json').read_text())
    assert info['feature_source']['engine'] == engine
    assert info['eigenvector_sign'] == 'solver-signs-then-final-vector-rule'


class MatlabFeatureTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name) / 'handoff'
        self.scans, self.metadata, self.regions = make_fixture()
        prepare(self.scans, self.metadata, self.regions, 2.1, self.folder)
        self.data = {
            'vectors': np.random.default_rng(35).normal(size=(135, 78)),
            'lengths': [40, 50, 45], 'tr': 2.1,
            'subjects': np.asarray([m['subject'] for m in self.metadata], dtype=object),
            'regions': np.asarray(self.regions, dtype=object),
            'input_signature': input_signature(self.scans, self.metadata, self.regions, 2.1),
            'runtime_engine': 'MATLAB', 'runtime_version': 'synthetic-test-only',
        }

    def test_prepare_uses_shared_preprocessing_once_and_preserves_order(self):
        exported = loadmat(self.folder / 'input.mat', simplify_cells=True)
        for actual, original in zip(exported['scans'], self.scans):
            np.testing.assert_array_equal(actual, preprocess_bold(original))
        self.assertEqual(exported['subjects'].tolist(), [m['subject'] for m in self.metadata])
        with self.assertRaises(FileExistsError):
            prepare(self.scans, self.metadata, self.regions, 2.1, self.folder)

    def test_import_runs_existing_pipeline_without_recomputing_features(self):
        savemat(self.folder / 'features.mat', self.data)
        check_export(self.folder, 'MATLAB')

    def test_mismatched_inputs_and_invalid_exports_are_rejected(self):
        path = self.folder / 'features.mat'
        corruptions = {'input_signature': 'wrong', 'subjects': self.data['subjects'][::-1],
                       'regions': self.data['regions'][::-1], 'lengths': [40, 49, 46], 'tr': 2.,
                       'vectors': np.full((135, 78), np.nan), 'runtime_engine': 'unknown'}
        for key, value in corruptions.items():
            with self.subTest(key=key):
                savemat(path, {**self.data, key: value})
                with self.assertRaises(ValueError):
                    load_features(path, self.scans, self.metadata, self.regions, 2.1)
        savemat(path, self.data)
        changed = [x.copy() for x in self.scans]
        changed[0][0, 0] += .01
        with self.assertRaises(ValueError):
            load_features(path, changed, self.metadata, self.regions, 2.1)


if __name__ == '__main__':
    unittest.main()
