"""共享后处理的数值、配对、缺失值与结果保存测试。"""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy import stats
from weida_postprocessing import FeatureTable, adjust_bh, compare_features, save_results
from weida_postprocessing.adapters import load_weida


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

    def test_bh_keeps_undefined_slots(self):
        np.testing.assert_allclose(adjust_bh([.02, .04, np.nan]), [.06, .06, np.nan])

    def test_shared_report_and_figures(self):
        with tempfile.TemporaryDirectory() as d:
            output = Path(d) / 'report'
            values = np.random.default_rng(4).normal(size=(9, 2))
            table = FeatureTable(values, self.scans, ('a', 'b'), ('one', 'one'), ('bit', 'bit'), 'test', 'test')
            save_results(table, output)
            self.assertEqual(len(list(output.glob('*.png'))), 2)
            self.assertEqual(len(list(output.glob('*.pdf'))), 2)
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
        self.assertEqual(result[0]['n_first'],2)
        self.assertEqual(result[0]['n_second'],2)
        self.assertEqual(result[1]['n_first'],3)
        for row in result: self.assertEqual(row['p'],row['q'])
        self.assertAlmostEqual(result[0]['p'],stats.ttest_rel(values[7:9,0],values[4:6,0]).pvalue)

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
