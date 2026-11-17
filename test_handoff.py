"""Check that handoff preparation is runtime-free and tampering is detected."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

class HandoffTests(unittest.TestCase):
    def test_prepare_and_checksum_guard(self):
        script=Path(__file__).parent/'validation'/'compare.py'
        with tempfile.TemporaryDirectory() as directory:
            result=subprocess.run([sys.executable,str(script),'--prepare-only','--output',directory],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            root=Path(directory)
            for name in ['input.mat','manifest.json','reference.m','run_reference.m','HANDOFF.md']:
                self.assertTrue((root/name).is_file(),name)
            with (root/'reference.m').open('a') as f: f.write('\n% changed\n')
            result=subprocess.run([sys.executable,str(script),'--compare-only','--output',directory],capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertIn('checksum mismatch',result.stderr)

if __name__=='__main__': unittest.main()
