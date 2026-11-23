"""Regression for numerical breakdown in the diagnostic candidate only."""
import unittest

import numpy as np

from validation.lanczos_probe import lanczos_probe
from weida import weight_eigenvectors


class LanczosProbeTests(unittest.TestCase):
    def test_rank_one_feature_has_known_value(self):
        for seed in (0, 13, 284):
            with self.subTest(seed=seed):
                start = np.random.default_rng(seed).normal(size=78)
                start /= np.linalg.norm(start)
                matrix = np.ones((78, 78))
                values, vectors, _ = lanczos_probe(matrix, start)
                np.testing.assert_allclose(
                    weight_eigenvectors(values, vectors),
                    -np.ones(78) / np.sqrt(78), atol=1e-12, rtol=0)

    def test_exact_invariant_subspace_stops_before_division_by_zero(self):
        matrix = np.diag(np.arange(1., 9.))
        start = np.eye(8)[:, -1]
        with np.errstate(all='raise'):
            values, vectors, details = lanczos_probe(matrix, start)
        self.assertTrue(details['breakdown'])
        self.assertEqual(details['steps'], 1)
        np.testing.assert_allclose(matrix @ vectors, vectors * values, atol=0)
