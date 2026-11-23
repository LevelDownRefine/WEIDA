"""Shared contracts: amplitude, edge orientation, study families and file roundtrip."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

from downstream import predict_states
from research_core.io import load_scans, read_saved_scan
from research_core.networks import module_summary
from research_core.preprocessing import preprocess_bold
from research_core.statistics import compare_feature_families


class SharedCoreTests(unittest.TestCase):
    def test_preprocessing_preserves_amplitude_and_rejects_invalid_input(self):
        time = np.arange(1000) * 2.1
        passband = np.sin(2 * np.pi * .055 * time)
        stopband = np.sin(2 * np.pi * .15 * time)
        scan = np.stack([passband, stopband])
        original = scan.copy()
        filtered = preprocess_bold(scan)
        np.testing.assert_array_equal(scan, original)
        np.testing.assert_allclose(preprocess_bold(3 * scan), 3 * filtered, atol=1e-10)
        self.assertGreater(np.std(filtered[0, 100:-100]), .6)
        self.assertLess(np.std(filtered[1, 100:-100]), .01)
        for x in (scan[:, :12], np.ones((2, 30)), np.full((2, 30), np.nan)):
            with self.assertRaises(ValueError):
                preprocess_bold(x)

    def test_module_summary_uses_both_directions_and_excludes_diagonal(self):
        matrix = np.array([[np.nan, 2, 1, 1], [6, np.nan, 1, 1],
                           [3, 3, np.nan, 4], [3, 3, 8, np.nan]])
        result = module_summary(matrix, [0, 0, 1, 1], directed=True)
        self.assertEqual(result['within'], 5)
        self.assertEqual(result['across'], 2)
        np.testing.assert_array_equal(result['nodal'], [1, 5, 1, 5])
        transposed = module_summary(matrix.T, [0, 0, 1, 1], directed=True)
        self.assertEqual(transposed['segregation'], result['segregation'])
        with self.assertRaisesRegex(ValueError, 'symmetric'):
            module_summary(matrix, [0, 0, 1, 1], directed=False)
        symmetric = (matrix + matrix.T) / 2
        symmetric_result = module_summary(symmetric, [0, 0, 1, 1], directed=False)
        self.assertEqual(symmetric_result['segregation'], 3)
        order = [2, 0, 3, 1]
        permuted = module_summary(matrix[np.ix_(order, order)], np.array([0, 0, 1, 1])[order], directed=True)
        self.assertEqual(permuted['segregation'], 3)
        np.testing.assert_array_equal(permuted['nodal'], result['nodal'][order])

    def test_feature_families_keep_subject_pairing_and_correction_separate(self):
        # One estimable feature, one undefined feature; preserve legacy finite-p BH.
        values = np.array([1, 2, 3, 2, 3, 5, 3, 5, 9.])
        features = {'scalar': values, 'with_missing': np.column_stack([values, np.full(9, np.nan)])}
        ids = ['a', 'b', 'c'] * 3
        groups = ['HC'] * 3 + ['OFF'] * 3 + ['ON'] * 3
        result = compare_feature_families(features, ids, groups)
        self.assertEqual(set(result), {'ON_vs_OFF', 'OFF_vs_HC', 'ON_vs_HC'})
        paired = result['ON_vs_OFF']['scalar']
        self.assertAlmostEqual(paired['mean_difference'][0], 7 / 3)
        self.assertAlmostEqual(paired['t'][0], np.sqrt(7))
        for contrast in result.values():
            self.assertEqual(contrast['scalar']['q'][0], contrast['with_missing']['q'][0])
            self.assertTrue(np.isnan(contrast['with_missing']['q'][1]))

    def test_downstream_cli_reloads_both_formats_with_identical_results(self):
        rng = np.random.default_rng(14)
        subjects = ['sub-HC01', 'sub-HC02', 'sub-HC03',
                    'sub-STNoff01', 'sub-STNoff02', 'sub-STNoff03',
                    'sub-STNon01', 'sub-STNon02', 'sub-STNon03']
        conditions = ['HC'] * 3 + ['STNoff'] * 3 + ['STNon'] * 3
        names = ['A', 'B', 'C', 'D']
        centers = np.array([[1, 1, -1, -1], [-1, -1, 1, 1.]])
        script = Path(__file__).resolve().parents[1] / 'run_features.py'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fmripost, pyleida = root / 'fmripost', root / 'pyleida'
            pyleida.mkdir()
            (pyleida / 'rois_labels.txt').write_text('\n'.join(names))
            with (pyleida / 'metadata.csv').open('w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['subject_id', 'condition'])
                writer.writerows(zip(subjects, conditions))
            for subject in subjects:
                x = rng.normal(size=(4, 30))
                folder = fmripost / subject
                folder.mkdir(parents=True)
                with (folder / (subject + '_BOLD_AAL1ctx.csv')).open('w', newline='') as handle:
                    writer = csv.writer(handle)
                    writer.writerow(['index'] + names)
                    writer.writerows([[i] + row.tolist() for i, row in enumerate(x.T)])
                folder = pyleida / 'time_series' / subject
                folder.mkdir(parents=True)
                np.savetxt(folder / 'tseries.csv', x, delimiter=',')
            vectors = rng.normal(size=(9 * 30, 4))
            lengths = np.full(9, 30)
            prediction = predict_states(vectors, centers, lengths)
            outputs = []
            for source in (fmripost, pyleida):
                scans, metadata, labels = load_scans(source, expected_regions=None)
                # Legacy model metadata must also reload correctly.
                legacy = {key: value for key, value in metadata[0].items() if key != 'scan_format'}
                np.testing.assert_array_equal(read_saved_scan(legacy, labels), scans[0])
                with self.assertRaisesRegex(ValueError, 'length changed'):
                    read_saved_scan(dict(legacy, frames=31), labels)
                model = root / (source.name + '_model')
                model.mkdir()
                (model / 'metadata.json').write_text(json.dumps({'scans': metadata, 'regions': labels, 'tr': 2.1}))
                np.savez(model / 'results.npz', vectors=vectors, lengths=lengths,
                         centers=centers, transition_valid=prediction['transitions'].sum(axis=2) > 0,
                         **prediction)
                output = root / (source.name + '_features')
                command = subprocess.run([sys.executable, str(script), str(model), '--output', str(output)],
                                         capture_output=True, text=True)
                self.assertEqual(command.returncode, 0, command.stderr)
                outputs.append(output)
            with np.load(outputs[0] / 'features.npz') as a, np.load(outputs[1] / 'features.npz') as b:
                for key in a.files:
                    # CSV loaders can return different memory layouts; reduction
                    # roundoff must not be mistaken for a scientific difference.
                    np.testing.assert_allclose(a[key], b[key], rtol=1e-12, atol=1e-14)
            for name in ('paired_statistics.json', 'independent_statistics.json'):
                self.assert_statistics_equal(json.loads((outputs[0] / name).read_text()),
                                             json.loads((outputs[1] / name).read_text()))

    def assert_statistics_equal(self, a, b):
        self.assertEqual(a.keys(), b.keys())
        for key in a:
            if isinstance(a[key], dict):
                self.assert_statistics_equal(a[key], b[key])
            elif key == 'paired_ids':
                self.assertEqual(a[key], b[key])
            else:
                np.testing.assert_allclose(a[key], b[key], rtol=1e-12, atol=1e-14)


if __name__ == '__main__':
    unittest.main()
