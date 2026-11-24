"""与算法无关的受试者指标表、统计检验及结果保存。"""
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests


@dataclass
class FeatureTable:
    values: np.ndarray
    scans: list
    names: tuple
    families: tuple
    units: tuple
    method: str
    definition: str

    def __post_init__(self):
        self.values=np.asarray(self.values,float)
        n=len(self.names)
        if self.values.shape!=(len(self.scans),n) or n==0 or np.isinf(self.values).any():
            raise ValueError('Expected scans x features, permitting NaN but not infinity')
        if len(set(self.names))!=n or len(self.families)!=n or len(self.units)!=n:
            raise ValueError('Feature names must be unique; every feature needs family and unit')
        if not self.method or not self.definition or not all(self.families):
            raise ValueError('Method, feature definition and correction families are required')
        group_indices(self.scans)


def group_indices(scans):
    groups={g:{} for g in ('HC','OFF','ON')}
    prefixes={'HC':'sub-HC','OFF':'sub-STNoff','ON':'sub-STNon'}
    for i,m in enumerate(scans):
        g=m['condition']
        if g not in groups: raise ValueError(f'Unsupported condition: {g}')
        if 'participant_id' in m:
            participant=str(m['participant_id'])
        else:
            if not m['subject'].startswith(prefixes[g]):
                raise ValueError('Supply participant_id when subject does not use DBS naming')
            participant=m['subject'][len(prefixes[g]):]
        if not participant or participant in groups[g]:
            raise ValueError('Empty or duplicate participant-condition')
        groups[g][participant]=i
    if any(len(v)<2 for v in groups.values()):
        raise ValueError('Require at least two participants per condition')
    if groups['OFF'].keys()!=groups['ON'].keys():
        raise ValueError('Incomplete ON/OFF pairs')
    return {g:[v[k] for k in sorted(v)] for g,v in groups.items()}


def adjust_bh(p):
    """原作者 get_fdr：直接调用 statsmodels BH，不替换 NaN。

    出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py
    """
    return multipletests(p, method="fdr_bh", alpha=0.05)[1]


def compare_features(table):
    """封装作者三组检验；仅增加 ID 对齐、指标族选择和具名输出。

    出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py
    不排除缺失值，不替换退化检验结果，不新增置信区间。
    """
    idx = group_indices(table.scans)
    rows = []
    for first, second in [('ON', 'OFF'), ('OFF', 'HC'), ('ON', 'HC')]:
        current = []
        paired = second == 'OFF'
        for k, name in enumerate(table.names):
            a, b = table.values[idx[first], k], table.values[idx[second], k]
            test = stats.ttest_rel(a, b) if paired else stats.ttest_ind(a, b)
            current.append({
                'contrast': first+'-'+second, 'feature': name, 'family': table.families[k],
                'unit': table.units[k], 'n_first': len(a), 'n_second': len(b),
                'mean_first': float(np.mean(a)), 'mean_second': float(np.mean(b)),
                'std_first': float(np.std(a)), 'std_second': float(np.std(b)),
                'difference': float(np.mean(a)-np.mean(b)),
                't': float(test.statistic), 'p': float(test.pvalue), 'q': np.nan,
                'test': 'paired_t' if paired else 'independent_equal_variance_t',
            })
        for family in dict.fromkeys(table.families):
            selected = [r for r in current if r['family'] == family]
            for row, q in zip(selected, adjust_bh([r['p'] for r in selected])):
                row['q'] = float(q)
        rows.extend(current)
    return rows


def write_table(path,rows):
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def save_results(table, output, *, input_path=None, features_path=None):
    """保存统计表、受试者指标表和元数据；防止覆盖已有分析结果。"""
    rows = compare_features(table)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    write_table(output / 'statistics.csv', rows)
    values = [dict(subject=scan['subject'], condition=scan['condition'],
                   **dict(zip(table.names, row))) for scan, row in zip(table.scans, table.values)]
    write_table(output / 'subject_features.csv', values)
    metadata = asdict(table)
    metadata.pop('values')
    metadata.update(input=str(Path(input_path).resolve()) if input_path is not None else None,
        features=str(Path(features_path).resolve()) if features_path is not None else None,
        correction='statsmodels fdr_bh within each contrast and declared family',
        missing='propagate; no deletion or replacement',
        reference='chenfei-ye/WEIDA-DBS@2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py')
    (output / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    return output
