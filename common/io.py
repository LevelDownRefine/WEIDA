"""PhiID 与 WEiDA 共用的 BOLD 读取、配对和表格保存。"""
import csv
from pathlib import Path

import numpy as np

__all__ = ['group_indices', 'load_scans', 'read_bold_csv', 'write_table']


def load_scans(root, exclude_stn=()):
    """加载现有 CSV 布局；标签与完整配对检查属于工程适配。"""
    scans, metadata, header = [], [], None
    root = Path(root)
    if (root / 'time_series').is_dir():
        header = (root / 'rois_labels.txt').read_text(encoding='utf-8-sig').splitlines()
        if not header or not all(header) or len(set(header)) != len(header):
            raise ValueError('Expected unique region labels')
        conditions = {'HC': 'HC', 'STNoff': 'OFF', 'STNon': 'ON'}
        with (root / 'metadata.csv').open(newline='', encoding='utf-8-sig') as handle:
            records = sorted(csv.DictReader(handle), key=lambda row: row['subject_id'])
        for row in records:
            subject = row['subject_id']
            condition = conditions[row['condition']]
            if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn):
                continue
            path = root / 'time_series' / subject / 'tseries.csv'
            x = np.loadtxt(path, delimiter=',', ndmin=2)
            if x.shape[0] != len(header):
                raise ValueError(f'Unexpected region count: {path}')
            scans.append(x)
            metadata.append({'subject': subject, 'condition': condition, 'file': str(path), 'frames': x.shape[1]})
    else:
        files = sorted(Path(root).glob('sub-*/*_BOLD_AAL1ctx*.csv'))
        for path in files:
            subject = path.parent.name
            if subject.startswith('sub-HC'): condition = 'HC'
            elif subject.startswith('sub-STNoff'): condition = 'OFF'
            elif subject.startswith('sub-STNon'): condition = 'ON'
            else: continue
            # 显式排除患者的 ON/OFF 两份扫描。
            if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn): continue
            x, names = read_bold_csv(path)
            if len(names) != 78 or (header is not None and names != header):
                raise ValueError(f'Unexpected atlas columns: {path}')
            header = names
            scans.append(x); metadata.append({'subject': subject, 'condition': condition, 'file': str(path), 'frames': x.shape[1]})
    if not scans: raise ValueError('No matching scans found')
    for x, record in zip(scans, metadata):
        if not x.size or not np.isfinite(x).all():
            raise ValueError(f"Empty or nonfinite scan: {record['file']}")
    off = {m['subject'].removeprefix('sub-STNoff') for m in metadata if m['condition']=='OFF'}
    on = {m['subject'].removeprefix('sub-STNon') for m in metadata if m['condition']=='ON'}
    if off != on: raise ValueError(f'Unpaired STN IDs: {off ^ on}')
    ids = [m['subject'] for m in metadata]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate scans for one subject/condition')
    return scans, metadata, header


def read_bold_csv(path):
    """作者 importdata → 删除首列 → 转置；保留全部时间点。

    出处：StateSpace_modelling.m 第 25–28、36–39 行。
    https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/StateSpace_modelling.m#L25-L39
    标签、形状及有限值检查为工程补充，防止脑区错位或坏数据静默进入分析。
    """
    path = Path(path)
    with path.open(newline='', encoding='utf-8-sig') as handle:
        rows = list(csv.reader(handle))
    if len(rows) < 2:
        raise ValueError(f'CSV 缺少表头或时间序列数据: {path}')
    names = rows[0][1:]
    if not names or not all(names) or len(set(names)) != len(names):
        raise ValueError(f'Expected unique atlas columns: {path}')
    if any(len(row) != len(names) + 1 for row in rows[1:]):
        raise ValueError(f'Unexpected CSV row width: {path}')
    values = np.array([[float(v) for v in row[1:]] for row in rows[1:]]).T
    if not np.isfinite(values).all():
        raise ValueError(f'Nonfinite scan: {path}')
    return values, names


def write_table(path, rows):
    """保存带字段名的结果表；PhiID 需记录分量、组别及数值，属于工程扩展。

    作者用 csvwrite 保存无表头数值矩阵（StateSpace_modelling.m 第 348–368 行）；
    本函数增加表头、文本列和 UTF-8 BOM，以保存 PhiID 的可解释结果。
    MATLAB 测试用 readtable 校验互操作，不宣称这是作者原函数。
    """
    rows = list(rows)
    if not rows:
        raise ValueError('Cannot save an empty result table')
    fields = list(rows[0])
    if not fields or any(set(row) != set(fields) for row in rows):
        raise ValueError('Result rows must have the same nonempty fields')
    with Path(path).open('w', newline='', encoding='utf-8-sig') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def group_indices(scans):
    """按 ID 配对；工程补充，替代作者固定人数及文件顺序假设。"""
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

