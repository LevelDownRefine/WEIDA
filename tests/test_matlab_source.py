"""Ensure the MATLAB reference remains byte-identical to the pinned upstream."""
import hashlib
import unittest
from pathlib import Path


class MatlabSourceTests(unittest.TestCase):
    def test_original_git_blobs(self):
        root = Path(__file__).parent / 'matlab' / 'upstream' / 'WEIDA-DBS'
        expected = {
            'StateSpace_modelling.m': '1d4b4a776d9c0bc99f509aa8d4867f6fefcae41c',
            'functions/demean.m': 'c5199f14864a19f9336b535d3667d604ebc8ee1f',
        }
        for name, digest in expected.items():
            with self.subTest(file=name):
                data = (root / name).read_bytes()
                blob = f'blob {len(data)}\0'.encode() + data
                self.assertEqual(hashlib.sha1(blob).hexdigest(), digest)
