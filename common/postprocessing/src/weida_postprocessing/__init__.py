"""与算法无关的受试者指标表、统计检验及结果保存。"""
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from . import author_statistics
from .author_statistics import (
    CONTRASTS,
    GROUPS,
    compare_groups,
    describe_groups,
    get_fdr,
    group_comparison,
)

__all__ = [
    "CONTRASTS",
    "GROUPS",
    "FeatureTable",
    "adjust_bh",
    "compare_features",
    "compare_groups",
    "describe_groups",
    "get_fdr",
    "group_comparison",
    "group_indices",
    "save_results",
    "write_table",
]


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


# 兼容已有名字，直接导出唯一的 BH 实现。
adjust_bh = author_statistics.get_fdr


def compare_features(table):
    """仅负责 ID 配对、指标族选择和输出；统计交给作者函数。"""
    idx = group_indices(table.scans)
    groups = tuple(table.values[idx[g]] for g in GROUPS)
    means, stds = author_statistics.describe_groups(*groups)
    t = np.empty((len(CONTRASTS), len(table.names)))
    p = np.empty_like(t)
    q = np.empty_like(t)
    for family in dict.fromkeys(table.families):
        columns = [k for k, value in enumerate(table.families) if value == family]
        t[:, columns], p[:, columns], q[:, columns] = author_statistics.group_comparison(
            *(values[:, columns] for values in groups))
    rows = []
    for contrast in ('ON-OFF', 'OFF-HC', 'ON-HC'):
        first, second = contrast.split('-')
        ci, ai, bi = CONTRASTS.index(contrast), GROUPS.index(first), GROUPS.index(second)
        for k, name in enumerate(table.names):
            rows.append({
                'contrast': contrast, 'feature': name, 'family': table.families[k],
                'unit': table.units[k], 'n_first': len(groups[ai]), 'n_second': len(groups[bi]),
                'mean_first': float(means[ai, k]), 'mean_second': float(means[bi, k]),
                'std_first': float(stds[ai, k]), 'std_second': float(stds[bi, k]),
                'difference': float(means[ai, k]-means[bi, k]),
                't': float(t[ci, k]), 'p': float(p[ci, k]), 'q': float(q[ci, k]),
                'test': 'paired_t' if second == 'OFF' else 'independent_equal_variance_t',
            })
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
