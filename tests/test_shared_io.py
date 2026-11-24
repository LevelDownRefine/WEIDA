"""公共读取、配对与表格输出；不依赖患者数据。"""
import csv
import tempfile
import unittest
from pathlib import Path

import numpy as np

from common.io import group_indices, read_bold_csv, write_table


class SharedIOTests(unittest.TestCase):
    def test_shape_validation_and_all_frames(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'bold.csv'
            path.write_text('index,A,B\n100,1,2\n101,3,4\n', encoding='utf-8-sig')
            values, names = read_bold_csv(path)
            np.testing.assert_array_equal(values, [[1, 3], [2, 4]])
            self.assertEqual(names, ['A', 'B'])
            for text in ['index,A,B\n0,1\n', 'index,A,A\n0,1,2\n', 'index,A\n0,nan\n']:
                path.write_text(text)
                with self.assertRaises(ValueError):
                    read_bold_csv(path)

    def test_pairing_reorders_ids_and_rejects_missing_pair(self):
        scans = [{'condition': g, 'participant_id': i} for g in ['HC', 'OFF', 'ON'] for i in ['2', '1']]
        indices = group_indices(scans)
        self.assertEqual(indices, {'HC': [1, 0], 'OFF': [3, 2], 'ON': [5, 4]})
        with self.assertRaises(ValueError):
            group_indices(scans[:-1])
        with self.assertRaises(ValueError):
            group_indices(scans + [scans[0]])

    def test_table_roundtrip_and_reject_before_writing(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'result.csv'
            rows = [{'atom': '冗余,"信息"', 'value': np.pi}, {'atom': 'rtr', 'value': -1e-20}]
            write_table(path, rows)
            with path.open(encoding='utf-8-sig', newline='') as handle:
                actual = list(csv.DictReader(handle))
            self.assertEqual(actual[0]['atom'], rows[0]['atom'])
            np.testing.assert_array_equal([float(row['value']) for row in actual], [np.pi, -1e-20])
            before = path.read_bytes()
            for invalid in [[], [{'a': 1}, {'b': 2}]]:
                with self.assertRaises(ValueError):
                    write_table(path, invalid)
                self.assertEqual(path.read_bytes(), before)
