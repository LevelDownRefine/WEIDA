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
        self.assertEqual(v.shape,(55,8))
        with self.assertRaises(ValueError): compute_weighted_eigenvectors([np.zeros((8,30))])
        c=cluster_brain_states(v,replicates=2)
        self.assertEqual(c['distances'].shape,(55,4))
        self.assertTrue(np.isfinite(c['silhouette_mean']))

if __name__=='__main__': unittest.main()
