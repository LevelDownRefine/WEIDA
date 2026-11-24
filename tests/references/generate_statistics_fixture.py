"""手动用固定版本作者原脚本生成基准；普通测试不联网。"""
import ast
import base64
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats
from statsmodels.stats.multitest import multipletests

REVISION = '2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c'


def main():
    result = subprocess.run(['gh', 'api', f'repos/chenfei-ye/WEIDA-DBS/contents/group_comparison.py?ref={REVISION}'],
                            check=True, capture_output=True, encoding='utf-8')
    response = json.loads(result.stdout)
    raw_source = base64.b64decode(response['content'])
    blob = hashlib.sha1(f'blob {len(raw_source)}\0'.encode() + raw_source).hexdigest()
    if blob != '19337d26053a985324f73b04bffb7e40262461e0':
        raise RuntimeError('原脚本内容与固定版本不一致')
    source = raw_source.decode('utf-8')
    node = next(n for n in ast.parse(source).body if isinstance(n, ast.FunctionDef) and n.name == 'get_fdr')
    fdr = ast.get_source_segment(source, node)
    block = source[source.index('df = np.zeros([4,12])'):]
    split = block.index('df_marker = ')
    inputs = np.random.default_rng(382).normal(size=(66, 4))
    cases = {}
    for name in ('finite', 'missing', 'constant', 'constant_change'):
        data = inputs.copy()
        if name == 'missing':
            data[16, 0] = np.nan
        elif name == 'constant':
            data[:, 0] = 1
        elif name == 'constant_change':
            data[16:41, 0] = 0
            data[41:, 0] = 1
        env = {'pd': pd, 'np': np, 'scipy': scipy, 'multipletests': multipletests,
               'data': pd.DataFrame(data), 'print': lambda *args: None}
        exec(fdr, env)  # noqa: S102 — 手动执行已校验的固定上游代码
        exec(block[:split], env)  # noqa: S102 — 原脚本计算段
        raw = env['df'].to_numpy(copy=True)
        exec(block[split:], env)  # noqa: S102 — 原脚本校正段
        corrected = env['df_marker'].to_numpy(copy=True)
        def serial(values):
            return [[float(v) if np.isfinite(v) else str(v) for v in row] for row in values]
        cases[name] = {key: serial(value) for key, value in {
            'mean': raw[:, [0, 2, 4]].T, 'std': raw[:, [1, 3, 5]].T,
            't': raw[:, [6, 8, 10]].T, 'p': raw[:, [7, 9, 11]].T,
            'q': corrected[:, [7, 9, 11]].T,
        }.items()}
    fixture = {'revision': REVISION, 'blob': response['sha'], 'input': inputs.tolist(), 'cases': cases}
    Path(__file__).with_name('statistics_fixture.json').write_text(
        json.dumps(fixture, ensure_ascii=False, indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__':
    main()
