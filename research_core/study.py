"""One public entry for loading, preprocessing, algorithm dispatch and inference."""
from dataclasses import dataclass
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np

from .backends import BACKENDS
from .contracts import AlgorithmResult
from .io import load_scans
from .preprocessing import preprocess_bold
from .statistics import compare_feature_families


def cohort_ids(metadata):
    groups = {key: set() for key in ('HC', 'OFF', 'ON')}
    prefixes = {'HC': 'sub-HC', 'OFF': 'sub-STNoff', 'ON': 'sub-STNon'}
    ids = []
    for row in metadata:
        condition = row['condition']
        if condition not in groups:
            raise ValueError(f'Unsupported condition: {condition}')
        participant = row.get('participant_id')
        if participant is None:
            if not row['subject'].startswith(prefixes[condition]):
                raise ValueError('Supply participant_id for non-DBS subject names')
            participant = row['subject'][len(prefixes[condition]):]
        participant = str(participant)
        if not participant or participant in groups[condition]:
            raise ValueError('Empty or duplicate participant-condition')
        groups[condition].add(participant)
        ids.append(participant)
    if groups['ON'] != groups['OFF']:
        raise ValueError('Incomplete ON/OFF pairs')
    if any(len(group) < 2 for group in groups.values()):
        raise ValueError('Require at least two participants per condition')
    return ids


@dataclass
class StudyResult:
    algorithm: AlgorithmResult
    scans: list
    regions: list
    comparisons: dict
    preprocessing: dict

    def save(self, output, *, plots=True):
        """Same filenames and metadata schema for every backend; never overwrite."""
        output = Path(output)
        output.mkdir(parents=True, exist_ok=False)
        np.savez_compressed(output / 'features.npz', **self.algorithm.features)
        np.savez_compressed(output / 'artifacts.npz', **self.algorithm.artifacts)
        metadata = dict(schema_version=1, method=self.algorithm.method, scans=self.scans,
                        regions=self.regions, preprocessing=self.preprocessing,
                        parameters=self.algorithm.parameters, names=self.algorithm.names,
                        units=self.algorithm.units, definitions=self.algorithm.definitions,
                        inference=dict(contrasts=['ON_vs_OFF', 'OFF_vs_HC', 'ON_vs_HC'],
                                       tests='two-sided paired / independent equal-variance t',
                                       correction='BH per contrast and family over finite p values'))
        for filename, data in [('metadata.json', metadata), ('statistics.json', self.comparisons)]:
            (output / filename).write_text(json.dumps(json_safe(data), ensure_ascii=False,
                                                     indent=2, allow_nan=False), encoding='utf-8')
        if plots:
            from .plotting import plot_study
            plot_study(self, output / 'figures')


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def analyze(scans, metadata, regions, *, method, preprocessing, tr=2.1, options=None):
    """Run either algorithm on the same region × time inputs.

    preprocessing is mandatory: 'loaded' or 'weida_amplitude'. Algorithm-specific
    options are strictly forwarded to that backend, never silently ignored.
    """
    if method not in BACKENDS:
        raise ValueError(f'Unknown method {method!r}; choose from {tuple(BACKENDS)}')
    if preprocessing not in ('loaded', 'weida_amplitude'):
        raise ValueError('Select loaded or weida_amplitude explicitly')
    if not np.isfinite(tr) or tr <= 0:
        raise ValueError('TR must be positive')
    scans = list(scans)
    metadata = [dict(row) for row in metadata]
    regions = list(regions)
    if len(scans) != len(metadata) or len(regions) < 2 or len(set(regions)) != len(regions):
        raise ValueError('Scan metadata or region labels do not match')
    ids = cohort_ids(metadata)
    options = dict(options or {})
    backend = BACKENDS[method]
    inspect.signature(backend).bind(scans, metadata, regions, tr=tr, **options)
    processed = []
    for scan, row in zip(scans, metadata):
        x = np.asarray(scan, dtype=float)
        if (x.ndim != 2 or x.shape[0] != len(regions) or x.shape[1] < 2
                or not np.isfinite(x).all() or np.any(np.ptp(x, axis=1) == 0)):
            raise ValueError(f"{row['subject']}: expected finite nonconstant region x time BOLD")
        row['frames'] = x.shape[1]
        row['input_array_sha256'] = hashlib.sha256(np.ascontiguousarray(x, dtype='<f8').tobytes()).hexdigest()
        processed.append(preprocess_bold(x, tr=tr) if preprocessing == 'weida_amplitude' else x.copy())
    algorithm = backend(processed, metadata, regions, tr=tr, **options).validate(len(scans))
    comparisons = compare_feature_families(algorithm.features, ids, [m['condition'] for m in metadata])
    profile = dict(name=preprocessing, tr=tr)
    if preprocessing == 'weida_amplitude':
        profile.update(detrend='linear', demean=True, band_hz=[.04, .07],
                       butterworth_order=2, filtfilt_padlen=12)
    return StudyResult(algorithm, metadata, regions, comparisons, profile)


def run_study(input_path, output, *, method, preprocessing, tr=2.1, exclude_stn=(),
              expected_regions=78, options=None, plots=True):
    """Filesystem entry; same loader, analysis, statistics and writer for both methods."""
    if Path(output).exists():
        raise FileExistsError(output)
    scans, metadata, regions = load_scans(input_path, exclude_stn, expected_regions=expected_regions)
    result = analyze(scans, metadata, regions, method=method, preprocessing=preprocessing,
                     tr=tr, options=options)
    result.preprocessing['excluded_stn'] = list(exclude_stn)
    result.save(output, plots=plots)
    return result
