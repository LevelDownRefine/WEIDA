"""薄入口：选择方法并调用现有脚本，不改变算法、统计或输出格式。"""
import argparse
import subprocess
import sys
from pathlib import Path


def run_study(input_path, output, *, method, phiid_repo=None, arguments=()):
    source, target = Path(input_path).resolve(), Path(output).resolve()
    if target.exists():
        raise FileExistsError(target)
    if method == 'weida':
        if phiid_repo is not None:
            raise ValueError('--phiid-repo only applies to PhiID')
        script = Path(__file__).resolve().with_name('run_weida.py')
        command = [sys.executable, str(script), str(source), str(target)]
    elif method == 'phiid':
        if phiid_repo is None:
            raise ValueError('Specify --phiid-repo pointing to the existing PhiID checkout')
        repo = Path(phiid_repo).resolve()
        script = repo / 'scripts' / 'run_study.py'
        if not script.is_file() or not (repo / 'uv.lock').is_file():
            raise ValueError(f'Expected PhiID checkout with scripts/run_study.py and uv.lock: {repo}')
        command = ['uv', 'run', '--locked', '--project', str(repo), 'python', str(script),
                   '--input', str(source), '--output', str(target)]
    else:
        raise ValueError(f'Unknown method: {method}')
    subprocess.run([*command, *arguments], check=True)
    return target


def main():
    parser = argparse.ArgumentParser(description=__doc__,
        epilog='在 -- 后传递原脚本参数，例如 -- --tr 2.1 --exclude-stn 18')
    parser.add_argument('--method', required=True, choices=['weida', 'phiid'])
    parser.add_argument('--phiid-repo', type=Path)
    parser.add_argument('input')
    parser.add_argument('output')
    # Explicit separator keeps method-specific flags under the original CLI's control.
    argv = sys.argv[1:]
    split = argv.index('--') if '--' in argv else len(argv)
    args = parser.parse_args(argv[:split])
    run_study(args.input, args.output, method=args.method, phiid_repo=args.phiid_repo,
              arguments=argv[split+1:])


if __name__ == '__main__':
    main()
