"""测试幅值预处理，以及进入方法专属相位提取前的边界。"""
import unittest
from unittest.mock import patch

import numpy as np
from weida_preprocessing import preprocess_scan

from preprocessing import preprocess_bold
from weida import compute_weighted_eigenvectors


class PreprocessingTests(unittest.TestCase):
    def test_package_preserves_opt_in_behavior(self):
        scan = np.arange(60.).reshape(2, 30)
        loaded = preprocess_scan(scan)
        np.testing.assert_array_equal(loaded, scan)
        self.assertFalse(np.shares_memory(loaded, scan))
        np.testing.assert_allclose(preprocess_scan(scan, demean=True).mean(axis=1), 0)
        np.testing.assert_allclose(preprocess_scan(scan, remove_trend=True), 0, atol=1e-12)

    def test_weida_wrapper_uses_package_result(self):
        scan = np.random.default_rng(284).normal(size=(78, 210))
        expected = preprocess_scan(scan, remove_trend=True, demean=True, band=(.04, .07), tr=2.1)
        np.testing.assert_array_equal(preprocess_bold(scan), expected)

    def test_frequency_response_amplitude_and_frame_count(self):
        time = np.arange(1000)*2.1
        scan = np.stack([np.sin(2*np.pi*.055*time), np.sin(2*np.pi*.15*time)])
        original = scan.copy()
        filtered = preprocess_bold(scan)
        self.assertEqual(filtered.shape, scan.shape)
        self.assertTrue(np.isrealobj(filtered))
        np.testing.assert_array_equal(scan, original)
        self.assertGreater(filtered[0, 100:-100].std(), .6)
        self.assertLess(filtered[1, 100:-100].std(), .01)
        np.testing.assert_allclose(preprocess_bold(3*scan), 3*filtered, atol=1e-10)

    def test_loaded_mode_copies_without_transforming(self):
        scan = np.arange(60.).reshape(2, 30)
        result = preprocess_bold(scan, band=None)
        np.testing.assert_array_equal(result, scan)
        self.assertFalse(np.shares_memory(result, scan))

    def test_staged_weida_matches_default_without_filtering_twice(self):
        rng = np.random.default_rng(37)
        scans = [rng.normal(size=(8, n)) for n in (80, 60)]
        expected, expected_lengths = compute_weighted_eigenvectors(scans)
        processed = [preprocess_bold(x) for x in scans]
        with patch('weida.preprocess_bold', side_effect=AssertionError('Repeated preprocessing')):
            actual, lengths = compute_weighted_eigenvectors(processed, preprocessed=True)
        np.testing.assert_array_equal(actual, expected)
        np.testing.assert_array_equal(lengths, expected_lengths)

    def test_invalid_inputs_fail_before_algorithm(self):
        scan = np.random.default_rng(15).normal(size=(3, 80))
        for kwargs in [{'tr': 0}, {'tr': np.nan}, {'band': (.07, .04)}, {'band': (.04, 1)}]:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                preprocess_bold(scan, **kwargs)
        for x in [scan[:, :12], np.ones((3, 30)), np.full((3, 30), np.nan), np.empty((3, 0))]:
            with self.assertRaises(ValueError):
                preprocess_bold(x)


if __name__ == '__main__':
    unittest.main()
