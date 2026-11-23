"""Thin launcher contract: delegation, argument boundaries and failure propagation."""
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from run_study import main, run_study


class LauncherTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)

    def test_weida_uses_existing_script_and_forwards_arguments(self):
        with patch('run_study.subprocess.run') as run:
            result = run_study(self.root/'input', self.root/'out', method='weida',
                               arguments=['--tr', '2.1', '--exclude-stn', '18'])
        command = run.call_args.args[0]
        self.assertEqual(command[0], sys.executable)
        self.assertEqual(Path(command[1]).name, 'run_weida.py')
        self.assertEqual(command[2:4], [str(self.root/'input'), str(self.root/'out')])
        self.assertEqual(command[4:], ['--tr', '2.1', '--exclude-stn', '18'])
        self.assertEqual(result, self.root/'out')
        self.assertEqual(run.call_args.kwargs, {'check': True})

    def test_phiid_uses_its_own_locked_environment(self):
        repo = self.root/'PhiID with spaces'
        (repo/'scripts').mkdir(parents=True)
        (repo/'scripts'/'run_study.py').touch()
        (repo/'uv.lock').touch()
        with patch('run_study.subprocess.run') as run:
            run_study(self.root/'input', self.root/'out', method='phiid', phiid_repo=repo,
                      arguments=['--preprocessing', 'loaded', '--tr', '2.1'])
        command = run.call_args.args[0]
        self.assertEqual(command[:7], ['uv', 'run', '--locked', '--project', str(repo),
                                      'python', str(repo/'scripts'/'run_study.py')])
        self.assertEqual(command[7:], ['--input', str(self.root/'input'), '--output', str(self.root/'out'),
                                      '--preprocessing', 'loaded', '--tr', '2.1'])

    def test_invalid_paths_or_methods_never_launch(self):
        with patch('run_study.subprocess.run') as run:
            for kwargs in [{'method': 'phiid'}, {'method': 'phiid', 'phiid_repo': self.root},
                           {'method': 'weida', 'phiid_repo': self.root}, {'method': 'unknown'}]:
                with self.assertRaises(ValueError):
                    run_study(self.root/'input', self.root/'out', **kwargs)
            with self.assertRaises(FileExistsError):
                run_study(self.root/'input', self.root, method='weida')
            run.assert_not_called()

    def test_child_failure_is_not_reported_as_success(self):
        with (patch('run_study.subprocess.run', side_effect=subprocess.CalledProcessError(2, 'child')),
              self.assertRaises(subprocess.CalledProcessError)):
            run_study(self.root/'input', self.root/'out', method='weida')

    def test_separator_preserves_backend_argument_boundaries(self):
        argv = ['run_study.py', '--method', 'phiid', '--phiid-repo', str(self.root),
                'input', 'output', '--', '--clinical', 'scores with spaces.csv', '--tr', '2.1']
        with patch.object(sys, 'argv', argv), patch('run_study.run_study') as run:
            main()
        self.assertEqual(run.call_args.kwargs['arguments'],
                         ['--clinical', 'scores with spaces.csv', '--tr', '2.1'])


if __name__ == '__main__':
    unittest.main()
