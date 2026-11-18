"""下游模块的手算与不变性测试。"""
import unittest
import numpy as np
from downstream import predict_states,state_metastability,functional_segregation,paired_comparison,independent_comparison

class DownstreamTests(unittest.TestCase):
    def test_prediction_preserves_center_order_and_tie_first(self):
        r=predict_states(np.array([[10.],[10.],[0.],[5.]]),np.array([[0.],[10.]]),[4])
        np.testing.assert_array_equal(r['states'],[1,1,0,0])
        np.testing.assert_allclose(r['transitions'][0],[[1,0],[.5,.5]])
        self.assertTrue(np.isnan(r['lifetime']).all())
    def test_independent_hand_t_and_missing_group(self):
        r=independent_comparison([0,2,2,4],['HC','HC','OFF','OFF'])
        self.assertAlmostEqual(r['mean_difference'][0],2)
        self.assertAlmostEqual(r['t'][0],np.sqrt(2))
        self.assertEqual(r['n_reference'][0],2)
        with self.assertRaises(ValueError): independent_comparison([0,1],['HC','HC'])
    def test_metastability_manual(self):
        r=state_metastability([[1,2],[3,4],[9,9]],np.array([0,0,1]),[[0,0],[9,9]],[3])
        self.assertEqual(r[0,0],15)
        self.assertTrue(np.isnan(r[0,1]))
    def test_segregation_orthogonal_signals(self):
        # 去均值且互相正交的四个向量，所有非对角相关为 0。
        x=np.array([[1,1,1,1,-1,-1,-1,-1],[1,1,-1,-1,1,1,-1,-1],[1,-1,1,-1,1,-1,1,-1],[1,-1,-1,1,-1,1,1,-1]],float)
        r=functional_segregation(x,[1,1,-1,-1])
        self.assertAlmostEqual(r['segregation'],0)
        np.testing.assert_allclose(r['nodal'],0,atol=1e-15)
        self.assertEqual(r['clipped_edges'],0)
    def test_perfect_correlations_are_counted(self):
        r=functional_segregation(np.tile(np.arange(10.),(4,1)),[1,1,-1,-1])
        self.assertEqual(r['clipped_edges'],6)
        self.assertTrue(np.isfinite(r['nodal']).all())
        with self.assertRaises(ValueError): functional_segregation(np.tile(np.arange(10.),(4,1)),[1,-1,-1,-1])
    def test_pairing_by_id_with_shuffled_rows_and_hand_t(self):
        r=paired_comparison([4,0,1,0,2,0],['c','a','a','c','b','b'],['ON','OFF','ON','OFF','ON','OFF'])
        self.assertEqual(r['paired_ids'],['a','b','c'])
        self.assertAlmostEqual(r['mean_difference'][0],7/3)
        self.assertAlmostEqual(r['t'][0],np.sqrt(7))
        self.assertEqual(r['n_pairs'][0],3)
        np.testing.assert_allclose(r['p'],r['q'])
    def test_missing_pair_and_duplicate_rejected(self):
        for ids,conditions in [(['a','b'],['OFF','ON']),(['a','a'],['OFF','OFF'])]:
            with self.assertRaises(ValueError): paired_comparison([0,1],ids,conditions)
    def test_nan_pair_exclusion_and_degenerate_difference(self):
        r=paired_comparison([[0,0],[1,1],[0,np.nan],[3,1]],['a','a','b','b'],['OFF','ON','OFF','ON'])
        np.testing.assert_array_equal(r['n_pairs'],[2,1])
        self.assertTrue(np.isnan(r['p'][1]))
        r=paired_comparison([0,1,0,1],['a','a','b','b'],['OFF','ON','OFF','ON'])
        self.assertTrue(np.isnan(r['p']).all())

if __name__=='__main__': unittest.main()
