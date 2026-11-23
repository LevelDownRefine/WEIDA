"""用解析解和手算预期值检验数值行为，而不是复制实现。"""
import unittest
from unittest.mock import patch

import numpy as np

from weida import (
    analyze_state_occupancy,
    calculate_state_dynamics,
    cluster_brain_states,
    compute_weighted_eigenvectors,
)


class NumericalOracleTests(unittest.TestCase):
    def test_identical_signals_have_known_rank_one_features(self):
        # 相同相位的同步矩阵全为 1；加权向量为负的单位常向量。
        wave=np.sin(np.arange(80)*.7)
        scans=np.tile(wave,(8,1)); original=scans.copy()
        features,lengths=compute_weighted_eigenvectors([scans])
        np.testing.assert_allclose(features,-np.ones((80,8))/np.sqrt(8),atol=1e-12)
        np.testing.assert_array_equal(scans,original)
        np.testing.assert_array_equal(lengths,[80])

    def test_preserve_individual_eigenvector_signs_before_weighting(self):
        # Inject solver output to isolate the original weighted-sum/sign rule.
        # Raw sum: [-1, 3]/4; only the final vector flips to [1, -3]/4.
        scan=np.tile(np.arange(4.),(2,1))
        with patch('weida.np.linalg.eigh',return_value=(np.array([1.,3.]),np.diag([-1.,1.]))):
            features,_=compute_weighted_eigenvectors([scan],preprocessed=True)
        np.testing.assert_array_equal(features,np.tile([.25,-.75],(4,1)))

    def test_filter_parameters_and_scan_compatibility(self):
        x=np.tile(np.arange(30),(8,1))
        for kwargs in [{'tr':0},{'tr':float('nan')},{'band':(.07,.04)},{'band':(.04,1)}]:
            with self.subTest(kwargs=kwargs),self.assertRaises(ValueError):
                compute_weighted_eigenvectors([x],**kwargs)
        for scans in [[],[x[:,:12]],[x,x[:7]]]:
            with self.assertRaises(ValueError): compute_weighted_eigenvectors(scans)

    def test_kmeans_distances_sumd_and_silhouette_by_hand(self):
        x=np.array([[0.],[2.],[10.],[12.]])
        result=cluster_brain_states(x,k=2,replicates=10)
        order=np.argsort(result['centers'][:,0])
        np.testing.assert_allclose(result['centers'][order],[[1.],[11.]])
        np.testing.assert_allclose(result['sumd'][order],[2.,2.])
        np.testing.assert_allclose(result['distances'][:,order],[[1,121],[1,81],[81,1],[121,1]])
        np.testing.assert_allclose(result['silhouette'],[118/122,78/82,78/82,118/122])

    def test_single_frame_has_no_observed_transitions(self):
        labels=np.array([0]); p,lt=analyze_state_occupancy(labels,[1],2)
        result=calculate_state_dynamics(labels,np.eye(2),[1],p,lt,['HC'])
        np.testing.assert_array_equal(result['transitions'],np.zeros((1,2,2)))
        self.assertTrue(np.isnan(result['group_transitions']['HC']).all())
        self.assertTrue(np.isnan(lt).all())

    def test_groups_do_not_mix_and_ties_keep_original_order(self):
        labels=np.array([0,0,1,1]); p,lt=analyze_state_occupancy(labels,[2,2],2)
        r=calculate_state_dynamics(labels,np.eye(2),[2,2],p,lt,['HC','ON'])
        np.testing.assert_array_equal(r['order'],[0,1])
        np.testing.assert_allclose(r['group_transitions']['HC'][0],[1,0])
        self.assertTrue(np.isnan(r['group_transitions']['HC'][1]).all())
        np.testing.assert_allclose(r['group_transitions']['ON'][1],[0,1])
        self.assertTrue(np.isnan(r['group_transitions']['ON'][0]).all())

if __name__=='__main__': unittest.main()
