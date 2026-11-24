"""直接验证已安装的正式统计实现与作者原脚本基准一致。"""
import json
import unittest
from pathlib import Path

import numpy as np
from weida_postprocessing import (
    CONTRASTS,
    GROUPS,
    compare_groups,
    describe_groups,
    group_comparison,
)


class AuthorStatisticsTests(unittest.TestCase):
    def test_refactor_preserves_original_script_results(self):
        fixture = json.loads((Path(__file__).parent / 'references/statistics_fixture.json').read_text(encoding='utf-8'))
        self.assertEqual(GROUPS, ('HC', 'OFF', 'ON'))
        self.assertEqual(CONTRASTS, ('OFF-HC', 'ON-OFF', 'ON-HC'))
        for case, expected in fixture['cases'].items():
            with self.subTest(case=case):
                data = np.array(fixture['input'])
                if case == 'missing':
                    data[16, 0] = np.nan
                elif case == 'constant':
                    data[:, 0] = 1
                elif case == 'constant_change':
                    data[16:41, 0] = 0
                    data[41:, 0] = 1
                groups = (data[:16], data[16:41], data[41:])
                means, stds = describe_groups(*groups)
                t, p, q = group_comparison(*groups)
                for key, actual in [('mean', means), ('std', stds), ('t', t), ('p', p), ('q', q)]:
                    np.testing.assert_allclose(actual, np.array(expected[key], dtype=float),
                                               rtol=0, atol=1e-12, equal_nan=True)
                raw_t, raw_p = compare_groups(*groups)
                np.testing.assert_allclose(raw_t, t, rtol=0, atol=0, equal_nan=True)
                np.testing.assert_allclose(raw_p, p, rtol=0, atol=0, equal_nan=True)
