"""共享后处理的数值、配对、缺失值与结果保存测试。"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import stats

from common.adapters import load_weida
from common.postprocessing import (
    FeatureTable,
    adjust_bh,
    compare_features,
    save_results,
    statistics,
)


class SharedResearchTests(unittest.TestCase):
    def setUp(self):
        self.scans=[{'subject': f'sub-{prefix}{i:02}','condition': g}
                    for g,prefix in [('HC','HC'),('OFF','STNoff'),('ON','STNon')] for i in range(1,4)]

    def test_same_statistics_regardless_of_method_name(self):
        values = np.random.default_rng(41).normal(size=(9, 2))
        def table(method):
            return FeatureTable(values, self.scans, ('a', 'b'), ('one', 'one'), ('u', 'u'), method, 'test')
        self.assertEqual(compare_features(table('PhiID')), compare_features(table('WEiDA')))
        result = compare_features(table('PhiID'))
        for k in range(2):
            self.assertAlmostEqual(result[k]['p'], stats.ttest_rel(values[6:, k], values[3:6, k]).pvalue)
            self.assertAlmostEqual(result[2+k]['p'], stats.ttest_ind(values[3:6, k], values[:3, k], equal_var=True).pvalue)
        order = [8, 2, 5, 0, 6, 3, 1, 7, 4]
        shuffled = FeatureTable(values[order], [self.scans[i] for i in order], ('a', 'b'),
                                ('one', 'one'), ('u', 'u'), 'test', 'test')
        self.assertEqual(compare_features(shuffled), result)

    def test_bh_preserves_author_nan_behavior(self):
        np.testing.assert_allclose(adjust_bh([.02, .04, np.nan]), np.full(3, np.nan), equal_nan=True)

    def test_shared_report_without_unvalidated_figures(self):
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / 'report'
            values = np.random.default_rng(4).normal(size=(9, 2))
            table = FeatureTable(values, self.scans, ('a', 'b'), ('one', 'one'), ('bit', 'bit'), 'test', 'test')
            save_results(table, output)
            self.assertEqual(len(list(output.glob('*.png'))), 0)
            self.assertEqual(len(list(output.glob('*.pdf'))), 0)
            self.assertTrue((output / 'statistics.csv').stat().st_size > 100)
            self.assertTrue((output / 'subject_features.csv').stat().st_size > 100)
            meta = json.loads((output / 'metadata.json').read_text(encoding='utf-8'))
            self.assertEqual(meta['units'], ['bit', 'bit'])
            self.assertEqual(meta['families'], ['one', 'one'])
            with self.assertRaises(FileExistsError):
                save_results(table, output)

    def test_missing_pairs_and_family_boundaries(self):
        values=np.random.default_rng(17).normal(size=(9,2));values[3,0]=np.nan
        table=FeatureTable(values,self.scans,('a','b'),('one','two'),('u','v'),'test','test')
        result=compare_features(table)
        self.assertEqual(result[0]['n_first'],3)
        self.assertEqual(result[0]['n_second'],3)
        self.assertEqual(result[1]['n_first'],3)
        for row in result: np.testing.assert_allclose(row['p'],row['q'],equal_nan=True)
        self.assertTrue(np.isnan(result[0]['p']))

    def test_adapter_matches_original_script_fixture(self):
        fixture = json.loads((Path(__file__).parent / 'references/statistics_fixture.json').read_text(encoding='utf-8'))
        scans = [{'subject': f'sub-{prefix}{i:02}', 'condition': group}
                 for group, prefix, n in [('HC', 'HC', 16), ('OFF', 'STNoff', 25), ('ON', 'STNon', 25)]
                 for i in range(n)]
        for case, expected in fixture['cases'].items():
            values = np.array(fixture['input'])
            if case == 'missing': values[16, 0] = np.nan
            if case == 'constant': values[:, 0] = 1
            if case == 'constant_change':
                values[16:41, 0] = 0
                values[41:, 0] = 1
            table = FeatureTable(values, scans, tuple('abcd'), ('phiid',)*4, ('bit',)*4, 'PhiID', 'test')
            rows = compare_features(table)
            for contrast, source in [('ON-OFF', 1), ('OFF-HC', 0), ('ON-HC', 2)]:
                selected = [r for r in rows if r['contrast'] == contrast]
                for key in ('t', 'p', 'q'):
                    np.testing.assert_allclose([r[key] for r in selected],
                        np.array(expected[key], dtype=float)[source], atol=1e-12, rtol=0, equal_nan=True)
                first, second = contrast.split('-')
                for suffix, group in [('first', first), ('second', second)]:
                    group_index = ('HC', 'OFF', 'ON').index(group)
                    for key in ('mean', 'std'):
                        np.testing.assert_allclose([r[f'{key}_{suffix}'] for r in selected],
                            np.array(expected[key], dtype=float)[group_index], atol=1e-12, rtol=0, equal_nan=True)
                self.assertNotIn('ci95_low', selected[0])

    def test_adapter_delegates_each_family_to_author_functions(self):
        from unittest.mock import patch
        table = FeatureTable(np.random.default_rng(15).normal(size=(9, 2)), self.scans,
                             ('a', 'b'), ('one', 'two'), ('u', 'u'), 'test', 'test')
        with patch.object(statistics, 'group_comparison', wraps=statistics.group_comparison) as test, \
             patch.object(statistics, 'describe_groups', wraps=statistics.describe_groups) as describe:
            compare_features(table)
            self.assertEqual(test.call_count, 2)
            self.assertEqual(describe.call_count, 1)
        self.assertIs(adjust_bh, statistics.get_fdr)

    def test_weida_saved_result_adapter(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);rng=np.random.default_rng(8)
            (root/'metadata.json').write_text(json.dumps({'scans': self.scans,'regions': ['a','b','c']}))
            valid=np.ones((9,2),bool);valid[0,1]=False
            np.savez(root/'results.npz',occupancy=rng.random((9,2)),lifetime=rng.random((9,2)),
                     transitions=rng.random((9,2,2)),transition_valid=valid)
            table=load_weida(root)
            self.assertEqual(table.values.shape,(9,8))
            self.assertTrue(np.isnan(table.values[0,-2:]).all())
            self.assertEqual(set(table.families),{'occupancy','lifetime','transitions'})
            self.assertEqual(table.units[2],'s')
            self.assertEqual(len(compare_features(table)),24)


if __name__=='__main__':unittest.main()
