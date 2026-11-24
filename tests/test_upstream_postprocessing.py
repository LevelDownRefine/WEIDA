"""与原作者 Python 后处理直接对照；正常输入等价与边界差异分别测试。"""
import unittest

import numpy as np
import pandas as pd

from downstream import functional_segregation, independent_comparison, paired_comparison
from tests.references.original_postprocessing import (
    calculate_functional_segregation_strength,
    calculate_nodal_functional_segregation_strength,
    get_fdr,
    original_group_comparison,
)


class UpstreamPostprocessingTests(unittest.TestCase):
    def test_three_group_statistics_and_bh(self):
        rng = np.random.default_rng(274)
        hc, off = rng.normal(size=(16, 4)), rng.normal(size=(25, 4))
        on = off + rng.normal(loc=.2, size=(25, 4))
        expected_t, expected_p, expected_q = original_group_comparison(hc, off, on)
        values = np.concatenate([hc, off, on])
        groups = ['HC'] * 16 + ['OFF'] * 25 + ['ON'] * 25
        ids = [f'h{i}' for i in range(16)] + list(range(25)) * 2
        actual = [independent_comparison(values, groups, target='OFF'),
                  paired_comparison(values, ids, groups),
                  independent_comparison(values, groups, target='ON')]
        for i, result in enumerate(actual):
            for key, expected in [('t', expected_t), ('p', expected_p), ('q', expected_q)]:
                np.testing.assert_allclose(result[key], expected[i], atol=1e-12, rtol=0)
        order = rng.permutation(len(values))
        shuffled = paired_comparison(values[order], [ids[i] for i in order], [groups[i] for i in order])
        for key in ('t', 'p', 'q'):
            np.testing.assert_allclose(shuffled[key], actual[1][key], atol=1e-12, rtol=0)

    def test_fdr_edge_values(self):
        from scipy.stats import false_discovery_control
        for p in ([0, 1, .01, .01], [.9, .9, .9, .9], [1e-12, .04, .2, .99]):
            np.testing.assert_allclose(false_discovery_control(p), get_fdr(p), atol=1e-12, rtol=0)

    def test_functional_and_nodal_segregation(self):
        rng = np.random.default_rng(275)
        for regions, frames, split in [(4, 80, 2), (7, 210, 2), (78, 210, 31)]:
            with self.subTest(regions=regions):
                scan = rng.normal(size=(regions, frames))
                center = np.r_[np.ones(split), -np.ones(regions-split)]
                ts = pd.DataFrame(scan.T)
                positive, negative = list(range(split)), list(range(split, regions))
                expected = calculate_functional_segregation_strength(ts, positive, negative)
                nodal = calculate_nodal_functional_segregation_strength(ts, positive, negative, regions)
                actual = functional_segregation(scan, center)
                np.testing.assert_allclose([actual[k] for k in ('within', 'across', 'segregation')], expected, atol=1e-12, rtol=0)
                np.testing.assert_allclose(actual['nodal'], nodal, atol=1e-12, rtol=0)
                self.assertEqual(actual['clipped_edges'], 0)

    def test_missing_values_are_an_explicit_extension(self):
        # 原脚本不排除 NaN，现实现按指标成对排除；此项记录差异，不宣称等价。
        off = np.array([[0.], [2.], [1.], [3.]])
        on = np.array([[1.], [4.], [np.nan], [7.]])
        _, p, _ = original_group_comparison(np.array([[0.], [1.], [2.]]), off, on)
        self.assertTrue(np.isnan(p[1, 0]))
        actual = paired_comparison(np.r_[off, on], list(range(4))*2, ['OFF']*4+['ON']*4)
        self.assertEqual(actual['n_pairs'][0], 3)
        self.assertTrue(np.isfinite(actual['p'][0]))


if __name__ == '__main__':
    unittest.main()
