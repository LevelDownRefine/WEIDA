import unittest

import numpy as np

from weida import *


class WorkflowTests(unittest.TestCase):
    def test_boundary_dwell(self):
        p,lt=analyze_state_occupancy(np.array([0,0,1,1,0,1]),[6],2,2)
        np.testing.assert_allclose(p,[[.5,.5]])
        np.testing.assert_allclose(lt,[[2,4]])
    def test_transitions_do_not_cross_scans(self):
        labels=np.array([0,0,1,1]); lengths=[2,2]
        p,lt=analyze_state_occupancy(labels,lengths,2)
        r=calculate_state_dynamics(labels,np.eye(2),lengths,p,lt,['A','A'])
        np.testing.assert_allclose(r['group_transitions']['A'],np.eye(2))
        np.testing.assert_allclose(r['transitions'][0,1],[0,0])
    def test_variable_lengths_and_rejection(self):
        x=np.random.default_rng(5).normal(size=(8,30))
        v,l=compute_weighted_eigenvectors([x,x[:,:25]])
        np.testing.assert_array_equal(l,[30,25])
        self.assertEqual(v.shape,(55,8))
        with self.assertRaises(ValueError): compute_weighted_eigenvectors([np.zeros((8,30))])
        c=cluster_brain_states(v,replicates=2)
        self.assertEqual(c['distances'].shape,(55,4))
        self.assertTrue(np.isfinite(c['silhouette_mean']))



class StateMetricsContractTests(unittest.TestCase):
    def test_absent_and_constant_states(self):
        p, lt = analyze_state_occupancy(np.array([1,1,1]), [3], 4)
        np.testing.assert_array_equal(p, [[0,1,0,0]])
        self.assertTrue(np.isnan(lt).all())

    def test_group_mean_is_scan_weighted_not_transition_pooled(self):
        labels=np.array([0,0,0,1, 0,1])
        p,lt=analyze_state_occupancy(labels,[4,2],2)
        r=calculate_state_dynamics(labels,np.eye(2),[4,2],p,lt,['A','A'])
        np.testing.assert_allclose(r['group_transitions']['A'][0],[1/3,2/3])
        self.assertTrue(np.isnan(r['group_transitions']['A'][1]).all())

    def test_state_reordering_is_consistent(self):
        labels=np.array([1,1,0,1])
        p,lt=analyze_state_occupancy(labels,[4],2)
        r=calculate_state_dynamics(labels,np.eye(2),[4],p,lt,['A'])
        np.testing.assert_array_equal(r['order'],[1,0])
        np.testing.assert_array_equal(r['states'],[0,0,1,0])
        np.testing.assert_allclose(r['occupancy'],[[.75,.25]])
        np.testing.assert_array_equal(r['centers'],np.eye(2)[[1,0]])

    def test_invalid_lengths_and_labels(self):
        for labels,lengths in [(np.array([0,1]),[3]),(np.array([0,4]),[2]),(np.array([0,1]),[0,2]),(np.array([0.,1.]),[2])]:
            with self.subTest(labels=labels,lengths=lengths), self.assertRaises(ValueError):
                analyze_state_occupancy(labels,lengths)

    def test_missing_and_nonfinite_regions_rejected(self):
        for x in [np.zeros((78,20)),np.full((78,20),np.nan),np.zeros((78,0))]:
            with self.assertRaises(ValueError): compute_weighted_eigenvectors([x])

if __name__=='__main__': unittest.main()
