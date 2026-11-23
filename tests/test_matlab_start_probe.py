"""Stream continuity tests; MATLAB sample comparison runs separately in CI."""
import unittest

import numpy as np

from validation.matlab_start_probe import DefaultStream, default_start


class MatlabStartProbeTests(unittest.TestCase):
    def test_chunking_does_not_change_the_normal_stream(self):
        stream = DefaultStream()
        chunked = np.concatenate([stream.normal(n) for n in (78, 1, 500, 1469)])
        np.testing.assert_array_equal(chunked, DefaultStream().normal(2048))

    def test_global_numpy_rng_does_not_change_default_start(self):
        state = np.random.get_state()
        try:
            start = default_start(78)
            np.random.seed(99)
            np.random.normal(size=100)
            np.testing.assert_array_equal(start, default_start(78))
            self.assertAlmostEqual(np.linalg.norm(start), 1.)
        finally:
            np.random.set_state(state)
