"""Verify shared pre/post processing and both real calculation methods."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from phyid.calculate import calc_PhiID

import run_study
from preprocessing import preprocess_bold
from weida import compute_weighted_eigenvectors


class StudyTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(84)
        self.scans = [rng.normal(size=(4, 48)) for _ in range(9)]
        self.meta = [{'subject': f'sub-{prefix}{i}', 'condition': g, 'frames': 48}
                     for prefix, g in [('HC', 'HC'), ('STNoff', 'OFF'), ('STNon', 'ON')] for i in range(3)]

    def test_only_method_changes_identical_preprocessed_input_and_postprocessing(self):
        captured, results = [], []
        def method(scans, conditions, tr):
            captured.append([x.copy() for x in scans])
            values = np.array([[x.mean(), x.std()] for x in scans])
            return {}, {'measurement': (values, ['mean', 'std'], 'a.u.')}, {}
        for name in ('weida', 'phiid'):
            with patch('run_study.'+name+'_method', side_effect=method):
                results.append(run_study.analyze(self.scans, self.meta, method=name))
        for a, b, raw in zip(captured[0], captured[1], self.scans):
            np.testing.assert_array_equal(a, b)
            np.testing.assert_array_equal(a, preprocess_bold(raw))
        for contrast in results[0][2]['measurement']:
            for key in ('t', 'p', 'q', 'mean_difference'):
                np.testing.assert_array_equal(results[0][2]['measurement'][contrast][key],
                                              results[1][2]['measurement'][contrast][key])

    def test_weida_preprocessing_once_matches_original_entry(self):
        before = [x.copy() for x in self.scans]
        with patch('run_study.preprocess_bold', wraps=preprocess_bold) as processing:
            native, _, _, _ = run_study.analyze(self.scans, self.meta, method='weida', k=2, replicates=2)
        self.assertEqual(processing.call_count, len(self.scans))
        expected, lengths = compute_weighted_eigenvectors(self.scans)
        np.testing.assert_array_equal(native['vectors'], expected)
        np.testing.assert_array_equal(native['lengths'], lengths)
        for a, b in zip(before, self.scans):
            np.testing.assert_array_equal(a, b)

    def test_official_phiid_numbers_and_region_order(self):
        native, features, _, _ = run_study.analyze(self.scans, self.meta, method='phiid', preprocessing='loaded')
        atoms, _ = calc_PhiID(self.scans[0][0], self.scans[0][1], 1, kind='gaussian', redundancy='MMI')
        expected = np.column_stack([atoms[name] for name in native['atom_names']]).mean(0)/np.log(2)
        np.testing.assert_allclose(native['atoms'][0, 0], expected, atol=1e-14)
        self.assertEqual(native['matrices'].shape, (9, 4, 4, 16))
        self.assertTrue(np.isnan(native['matrices'][:, np.arange(4), np.arange(4)]).all())
        _, permuted, _, _ = run_study.analyze([x[[2, 0, 3, 1]] for x in self.scans], self.meta,
                                            method='phiid', preprocessing='loaded')
        np.testing.assert_allclose(features['phiid_mean'][0], permuted['phiid_mean'][0], atol=1e-13)
        names = list(native['atom_names'])
        swap = [names.index(a.translate(str.maketrans('xy', 'yx'))) for a in names]
        np.testing.assert_allclose(features['phiid_mean'][0], features['phiid_mean'][0][:, swap], atol=1e-14)

    def test_same_loader_writer_and_plots_for_both_real_methods(self):
        with tempfile.TemporaryDirectory() as temporary:
            for method, options in [('weida', {'k': 2, 'replicates': 2}), ('phiid', {'lag': 1})]:
                output = Path(temporary)/method
                with patch('run_study.load_scans', return_value=(self.scans, self.meta, ['A', 'B', 'C', 'D'])) as loader:
                    run_study.run_study('input', output, method=method, exclude_stn=['18'], **options)
                    loader.assert_called_once_with('input', ['18'])
                info = json.loads((output/'metadata.json').read_text(encoding='utf-8'))
                self.assertEqual(info['preprocessing'], 'weida_amplitude')
                self.assertEqual(info['method'], method)
                self.assertEqual(info['excluded_stn'], ['18'])
                self.assertEqual(len(list(output.glob('*.png'))), 2*len(info['features']))
                self.assertEqual(len(list(output.glob('*.pdf'))), 2*len(info['features']))
                self.assertNotIn('NaN', (output/'statistics.json').read_text())
                with self.assertRaises(FileExistsError):
                    run_study.run_study('input', output, method=method)

    def test_invalid_input_and_inapplicable_parameters_fail(self):
        for kwargs in [{'method': 'unknown'}, {'method': 'weida', 'preprocessing': 'phase'},
                       {'method': 'phiid', 'tr': 0}, {'method': 'phiid', 'lag': 0}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                run_study.analyze(self.scans, self.meta, **kwargs)
        with self.assertRaises(TypeError):
            run_study.analyze(self.scans, self.meta, method='phiid', k=2)
        with self.assertRaisesRegex(ValueError, 'complete'):
            run_study.analyze(self.scans[:-1], self.meta[:-1], method='phiid')
        singular = [x.copy() for x in self.scans]
        singular[0][1] = singular[0][0]
        with self.assertRaisesRegex(ValueError, 'singular'):
            run_study.analyze(singular, self.meta, method='phiid', preprocessing='loaded')


if __name__ == '__main__':
    unittest.main()
