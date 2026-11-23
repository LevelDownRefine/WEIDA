"""Read the PD DBS cohort without algorithm-specific transformations."""
import csv
from pathlib import Path
import numpy as np


def read_scan(path, regions, *, scan_format):
    """Read a region × time array, checking labels for header-bearing CSVs."""
    path = Path(path)
    if scan_format == 'pyleida':
        x = np.loadtxt(path, delimiter=',', ndmin=2)
    elif scan_format == 'fmripost':
        with path.open(newline='', encoding='utf-8-sig') as handle:
            rows = list(csv.reader(handle))
        if len(rows) < 2:
            raise ValueError(f'CSV 缺少表头或时间序列数据: {path}')
        if rows[0][1:] != list(regions):
            raise ValueError(f'Unexpected atlas columns: {path}')
        x = np.array([[float(v) for v in row[1:]] for row in rows[1:]]).T
    else:
        raise ValueError(f'Unknown scan format: {scan_format}')
    if x.ndim != 2 or x.shape[0] != len(regions) or x.shape[1] == 0 or not np.isfinite(x).all():
        raise ValueError(f'Expected nonempty finite region x time array: {path}')
    return x


def read_saved_scan(metadata, regions):
    """Reload a scan, including outputs written before scan_format was recorded.

    Legacy PyLEiDA files are identified only by their known tseries.csv name.
    Headerless files rely on the dataset-wide region order saved in metadata.
    """
    path = Path(metadata['file'])
    scan_format = metadata.get('scan_format', 'pyleida' if path.name == 'tseries.csv' else 'fmripost')
    x = read_scan(path, regions, scan_format=scan_format)
    if x.shape[1] != metadata['frames']:
        raise ValueError(f'Scan length changed since model fitting: {path}')
    return x


def load_scans(root, exclude_stn=(), *, expected_regions=78):
    """Load HC/STN scans with strict pairing; expected_regions=None accepts other atlases."""
    scans, metadata, header = [], [], None
    root = Path(root)
    if (root / 'time_series').is_dir():
        header = (root / 'rois_labels.txt').read_text(encoding='utf-8-sig').splitlines()
        if (not header or len(set(header)) != len(header) or not all(header)
                or (expected_regions is not None and len(header) != expected_regions)):
            raise ValueError('Unexpected atlas labels or region count')
        conditions = {'HC': 'HC', 'STNoff': 'OFF', 'STNon': 'ON'}
        with (root / 'metadata.csv').open(newline='', encoding='utf-8-sig') as handle:
            records = sorted(csv.DictReader(handle), key=lambda row: row['subject_id'])
        for row in records:
            subject = row['subject_id']
            condition = conditions[row['condition']]
            if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn):
                continue
            path = root / 'time_series' / subject / 'tseries.csv'
            x = read_scan(path, header, scan_format='pyleida')
            scans.append(x)
            metadata.append({'subject': subject, 'condition': condition, 'file': str(path.resolve()), 'frames': x.shape[1], 'scan_format': 'pyleida'})
    else:
        files = sorted(Path(root).glob('sub-*/*_BOLD_AAL1ctx*.csv'))
        for path in files:
            subject = path.parent.name
            if subject.startswith('sub-HC'): condition = 'HC'
            elif subject.startswith('sub-STNoff'): condition = 'OFF'
            elif subject.startswith('sub-STNon'): condition = 'ON'
            else: continue
            # Exclude the invalid scan and its pair explicitly.
            if any(subject in ('sub-STNon'+n, 'sub-STNoff'+n) for n in exclude_stn): continue
            with path.open(newline='', encoding='utf-8-sig') as f:
                rows = list(csv.reader(f))
            if len(rows) < 2:
                raise ValueError(f'CSV 缺少表头或时间序列数据: {path}')
            names = rows[0][1:]
            if (not names or not all(names) or len(set(names)) != len(names)
                    or (expected_regions is not None and len(names) != expected_regions)
                    or (header is not None and names != header)):
                raise ValueError(f'Unexpected atlas columns: {path}')
            header = names
            x = read_scan(path, header, scan_format='fmripost')
            scans.append(x); metadata.append({'subject': subject, 'condition': condition, 'file': str(path.resolve()), 'frames': x.shape[1], 'scan_format': 'fmripost'})
    if not scans: raise ValueError('No matching scans found')
    off = {m['subject'].removeprefix('sub-STNoff') for m in metadata if m['condition']=='OFF'}
    on = {m['subject'].removeprefix('sub-STNon') for m in metadata if m['condition']=='ON'}
    if off != on: raise ValueError(f'Unpaired STN IDs: {off ^ on}')
    ids = [m['subject'] for m in metadata]
    if len(set(ids)) != len(ids): raise ValueError('Duplicate scans for one subject/condition')
    return scans, metadata, header
