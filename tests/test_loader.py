"""CSV 入口测试：临时合成文件，不依赖服务器数据。"""
import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
from run_weida import load_scans

class LoaderTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root=Path(self.temp.name)
        self.names=[f'ROI{i}' for i in range(78)]

    def write_scan(self,subject,suffix='',names=None):
        folder=self.root/subject; folder.mkdir(exist_ok=True)
        path=folder/f'{subject}_BOLD_AAL1ctx{suffix}.csv'
        values=np.arange(234).reshape(3,78)
        with path.open('w',newline='',encoding='utf-8-sig') as f:
            writer=csv.writer(f); writer.writerow(['']+(self.names if names is None else names))
            for i,row in enumerate(values): writer.writerow([100+i]+row.tolist())
        return path,values

    def test_bom_index_removed_and_time_axis_transposed(self):
        _,values=self.write_scan('sub-HC01')
        scans,meta,names=load_scans(self.root)
        np.testing.assert_array_equal(scans[0],values.T)
        self.assertEqual(meta[0]['frames'],3)
        self.assertEqual(meta[0]['condition'],'HC')
        self.assertEqual(names,self.names)

    def test_pyleida_input_matches_fmripost_including_first_frame(self):
        _, values = self.write_scan('sub-HC01')
        expected, _, labels = load_scans(self.root)
        root = self.root / 'input'
        folder = root / 'time_series' / 'sub-HC01'
        folder.mkdir(parents=True)
        (root / 'rois_labels.txt').write_text('\n'.join(labels))
        (root / 'metadata.csv').write_text('subject_id,condition\nsub-HC01,HC\n')
        np.savetxt(folder / 'tseries.csv', values.T, delimiter=',')
        scans, meta, names = load_scans(root)
        np.testing.assert_array_equal(scans[0], expected[0])
        self.assertEqual(names, labels)
        self.assertEqual(meta[0]['frames'], 3)

    def test_explicit_exclusion_removes_both_scans_and_ignores_gpi(self):
        for subject in ['sub-HC01','sub-STNon18','sub-STNoff18','sub-GPIon01']:
            self.write_scan(subject)
        _,meta,_=load_scans(self.root,exclude_stn=['18'])
        self.assertEqual([m['subject'] for m in meta],['sub-HC01'])

    def test_pairs_match_ids_not_only_counts(self):
        self.write_scan('sub-STNon01'); self.write_scan('sub-STNoff02')
        with self.assertRaisesRegex(ValueError,'Unpaired'): load_scans(self.root)

    def test_duplicate_scan_is_rejected(self):
        self.write_scan('sub-HC01','_a'); self.write_scan('sub-HC01','_b')
        with self.assertRaisesRegex(ValueError,'Duplicate'): load_scans(self.root)

    def test_region_order_mismatch_is_rejected(self):
        self.write_scan('sub-HC01'); self.write_scan('sub-HC02',names=self.names[::-1])
        with self.assertRaisesRegex(ValueError,'atlas'): load_scans(self.root)

    def test_empty_directory_is_rejected(self):
        with self.assertRaisesRegex(ValueError,'No matching'): load_scans(self.root)

    def test_empty_csv_has_contextual_error(self):
        path,_=self.write_scan('sub-HC01'); path.write_text('')
        with self.assertRaisesRegex(ValueError,'sub-HC01'): load_scans(self.root)

    def test_duplicate_region_names_are_rejected(self):
        self.write_scan('sub-HC01',names=['same']*78)
        with self.assertRaises(ValueError): load_scans(self.root)

if __name__=='__main__': unittest.main()
