"""Execute both real backends through the same public API and CLI."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from phyid.calculate import calc_PhiID

from research_core.backends import ATOM_NAMES, SWAP_XY
from research_core.preprocessing import preprocess_bold
from research_core.study import analyze, run_study
from weida import (compute_weighted_eigenvectors, cluster_brain_states,
                   analyze_state_occupancy, calculate_state_dynamics)


class UnifiedStudyTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(81)
        cls.scans = [rng.normal(size=(4, 64)) for _ in range(9)]
        cls.metadata = [dict(subject=prefix+str(i), condition=condition)
                        for prefix, condition in [('sub-HC', 'HC'), ('sub-STNoff', 'OFF'), ('sub-STNon', 'ON')]
                        for i in range(1, 4)]
        cls.regions = ['A', 'B', 'C', 'D']

    def test_weida_shared_preprocessing_runs_once_and_matches_legacy(self):
        originals = [x.copy() for x in self.scans]
        with patch('research_core.study.preprocess_bold', wraps=preprocess_bold) as preprocessing:
            result = analyze(self.scans, self.metadata, self.regions, method='weida',
                             preprocessing='weida_amplitude', options=dict(k=2, seed=3, replicates=3))
            self.assertEqual(preprocessing.call_count, len(self.scans))
        vectors, lengths = compute_weighted_eigenvectors(self.scans)
        clustering = cluster_brain_states(vectors, k=2, random_state=3, replicates=3)
        p, lt = analyze_state_occupancy(clustering['labels'], lengths, k=2)
        expected = calculate_state_dynamics(clustering['labels'], clustering['centers'], lengths,
                                            p, lt, [m['condition'] for m in self.metadata])
        np.testing.assert_array_equal(result.algorithm.artifacts['vectors'], vectors)
        for key in ('states', 'centers', 'occupancy', 'lifetime', 'transitions', 'transition_valid'):
            if key == 'centers':
                # Parallel K-means centroid reductions may differ at rounding precision.
                np.testing.assert_allclose(result.algorithm.artifacts[key], expected[key], atol=1e-14, rtol=1e-12)
            else:
                np.testing.assert_array_equal(result.algorithm.artifacts[key], expected[key])
        for actual, original in zip(self.scans, originals):
            np.testing.assert_array_equal(actual, original)

    def test_phiid_calls_official_estimator_and_records_units_and_lag(self):
        result = analyze(self.scans, self.metadata, self.regions, method='phiid',
                         preprocessing='loaded', options={'lag': 2})
        atoms, _ = calc_PhiID(self.scans[0][0], self.scans[0][1], 2, kind='gaussian', redundancy='MMI')
        expected = np.column_stack([atoms[name] for name in ATOM_NAMES]).mean(0) / np.log(2)
        np.testing.assert_allclose(result.algorithm.artifacts['atoms'][0, 0], expected, atol=1e-14)
        matrices = result.algorithm.artifacts['matrices']
        self.assertEqual(matrices.shape, (9, 4, 4, 16))
        np.testing.assert_allclose(matrices[0, 1, 0], expected[SWAP_XY], atol=1e-14)
        self.assertTrue(np.isnan(matrices[:, np.arange(4), np.arange(4)]).all())
        self.assertEqual(result.algorithm.units['phiid_mean'], 'bit')
        self.assertEqual(result.algorithm.parameters['lag_seconds'], 4.2)
        self.assertEqual(set(result.comparisons), {'ON_vs_OFF', 'OFF_vs_HC', 'ON_vs_HC'})

    def test_phiid_mean_invariant_to_region_order_and_xy_swap(self):
        first = analyze(self.scans, self.metadata, self.regions, method='phiid', preprocessing='loaded')
        order = [2, 0, 3, 1]
        second = analyze([x[order] for x in self.scans], self.metadata,
                         [self.regions[i] for i in order], method='phiid', preprocessing='loaded')
        values = first.algorithm.features['phiid_mean']
        np.testing.assert_allclose(values, values[:, SWAP_XY], atol=1e-14)
        np.testing.assert_allclose(values, second.algorithm.features['phiid_mean'], atol=1e-13)
        for contrast in first.comparisons:
            np.testing.assert_allclose(first.comparisons[contrast]['phiid_mean']['p'],
                                       second.comparisons[contrast]['phiid_mean']['p'], atol=1e-10)

    def test_invalid_dispatch_cohort_and_backend_options_fail(self):
        for kwargs in [dict(method='unknown', preprocessing='loaded'),
                       dict(method='phiid', preprocessing='phase'),
                       dict(method='phiid', preprocessing='loaded', tr=0),
                       dict(method='phiid', preprocessing='loaded', options={'lag': 0})]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                analyze(self.scans, self.metadata, self.regions, **kwargs)
        with self.assertRaises(TypeError):
            analyze(self.scans, self.metadata, self.regions, method='phiid',
                    preprocessing='loaded', options={'k': 2})
        with self.assertRaisesRegex(ValueError, 'Incomplete'):
            analyze(self.scans[:-1], self.metadata[:-1], self.regions, method='weida', preprocessing='loaded')
        singular = [x.copy() for x in self.scans]
        singular[0][1] = singular[0][0]
        with self.assertRaisesRegex(ValueError, 'sub-HC1.*A/B'):
            analyze(singular, self.metadata, self.regions, method='phiid', preprocessing='loaded')

    def test_same_cli_dispatches_both_methods_and_exports_common_schema(self):
        script = Path(__file__).resolve().parents[1] / 'run_study.py'
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'rois_labels.txt').write_text('\n'.join(self.regions))
            with (root / 'metadata.csv').open('w', newline='') as handle:
                writer = csv.writer(handle)
                writer.writerow(['subject_id', 'condition'])
                for row, scan in zip(self.metadata, self.scans):
                    writer.writerow([row['subject'], {'HC': 'HC', 'OFF': 'STNoff', 'ON': 'STNon'}[row['condition']]])
                    folder = root / 'time_series' / row['subject']
                    folder.mkdir(parents=True)
                    np.savetxt(folder / 'tseries.csv', scan, delimiter=',')
            for method, options in [('weida', ['--k', '2', '--replicates', '2']), ('phiid', ['--lag', '1'])]:
                output = root / method
                command = subprocess.run([sys.executable, str(script), str(root), str(output),
                                          '--method', method, '--preprocessing', 'weida_amplitude',
                                          '--expected-regions', '4', '--no-plots', *options],
                                         capture_output=True, text=True)
                self.assertEqual(command.returncode, 0, command.stderr)
                metadata = json.loads((output / 'metadata.json').read_text(encoding='utf-8'))
                self.assertEqual(metadata['method'], method)
                self.assertEqual(metadata['preprocessing']['name'], 'weida_amplitude')
                self.assertEqual(len(metadata['scans']), 9)
                self.assertIn('BH per contrast', metadata['inference']['correction'])
                for filename in ('features.npz', 'artifacts.npz', 'statistics.json'):
                    self.assertTrue((output / filename).is_file())
                self.assertNotIn('NaN', (output / 'statistics.json').read_text())
                with self.assertRaises(FileExistsError):
                    run_study(root, output, method=method, preprocessing='loaded')

    def test_shared_writer_and_plotter_accept_both_backends(self):
        with tempfile.TemporaryDirectory() as directory:
            for method, options in [('weida', dict(k=2, replicates=2)), ('phiid', dict(lag=1))]:
                result = analyze(self.scans, self.metadata, self.regions, method=method,
                                 preprocessing='loaded', options=options)
                output = Path(directory) / method
                result.save(output)
                pngs = list((output / 'figures').glob('*.png'))
                self.assertEqual(len(pngs), 2 * len(result.algorithm.features))
                for png in pngs:
                    self.assertGreater(png.stat().st_size, 1000)
                    self.assertTrue(png.with_suffix('.pdf').is_file())


if __name__ == '__main__':
    unittest.main()
