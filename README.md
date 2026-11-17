# WEiDA Python

Four functions mirror the four numbered stages of [StateSpace_modelling.m](https://github.com/chenfei-ye/WEIDA-DBS/blob/main/StateSpace_modelling.m):

1. `compute_weighted_eigenvectors`: detrend, demean, 0.04–0.07 Hz Butterworth filtering, Hilbert phase, weighted eigenvectors.
2. `cluster_brain_states`: four-state K-means (300 restarts, 500 iterations), squared distances, within-cluster sums and exact squared-Euclidean silhouette.
3. `analyze_state_occupancy`: occupancy and dwell time excluding incomplete boundary runs.
4. `calculate_state_dynamics`: order states by overall frequency, calculate scan transitions and eligible-scan group averages.

## Run

```
uv sync --locked
uv run --locked python run_weida.py /path/to/PD_DBS_zcc/derivatives/fmripost output --exclude-stn 18 --tr 2.1
uv run --locked python -m unittest discover -v
```

The current dataset has 15 HC and 26 STN pairs. STNon18 has 53 zero regions; explicitly excluding pair 18 gives 15 HC + 25 OFF + 25 ON scans. The loader verifies matching atlas headers, unique scans and ON/OFF pairing. GPI is not included. Exclusions are recorded; no subjects are silently excluded. Verify actual TR and scientific inclusion criteria before interpreting results. The input is never modified.

## Numerical contract and limitations

This is an equation-level migration, NOT a validated MATLAB reproduction. The source weights multiple eigenvectors before fixing the sign of their sum. Individual eigenvectors have arbitrary signs, so swapping eigensolvers can change the scientific result. This implementation explicitly orients EACH eigenvector so its largest absolute component is positive, then applies the source's overall sign rule. This is an intentional methodological convention, not a demonstrated equivalence to MATLAB eigs. Near-degenerate eigenspaces remain ambiguous. A MATLAB reference comparison is still required before scientific use.

Filtering explicitly uses padding length `3*(max(len(a),len(b))-1)` rather than SciPy's different default; MATLAB version-specific endpoint behavior remains unverified. K-means random initialization differs across libraries. Results must be compared after state-label matching, not by raw label numbers. The original silhouette is calculated exactly and may be expensive for large datasets.

## Outputs

`clustering.npz` retains unsorted labels, centers, squared distances, SUMD-equivalent sums and silhouette. `results.npz` contains sorted centers, occupancy, dwell time, zero-based states, scan transitions, valid-source masks and group transitions. Scan boundaries are stored in `lengths`; `metadata.json` records scan order, region names and parameters.

CSV outputs include `Vemp.csv` (regions x states), `State_ts.csv` (one-based), `Ptrans1.csv` etc., and HC/OFF/ON group files `P1emp.csv`, `LT1emp.csv`, `PTR1emp.csv` etc. Occupancy/dwell CSVs have only the corresponding group's rows, unlike MATLAB's all-scan NaN placeholders; use metadata to recover subject IDs. Missing source states yield zero rows per scan, as in MATLAB; group averaging excludes those scans for that source state. No valid contributors yields NaN.

Tests cover complete-run dwell, scan-boundary isolation, transition averaging, unequal lengths, zero-series rejection and clustering outputs. Full real-data execution and MATLAB reference validation have not yet been performed. Upstream code and paper are credited above; no upstream license or permission to redistribute original source is asserted here.

## Portable Octave validation

See [validation/README.md](validation/README.md) for a runnable synthetic cross-language comparison and its recorded report. Octave 11.3.0 passed 16 of 17 comparisons; the original-sign weighted eigenvectors differ, as expected from the explicitly changed sign convention. This does not establish MATLAB equivalence.

## Environment

The project uses uv with Python 3.11.14 (`.python-version`). Dependencies are declared in `pyproject.toml` and resolved exactly in `uv.lock`; `.venv` stays local. Run all commands above from the repository root. CI uses the same Python pin and `uv sync --locked`, so it fails if dependency declarations and the lockfile disagree. Update dependencies deliberately with `uv lock --upgrade`, then rerun unit tests and the Octave comparison before committing the new lockfile.

MATLAB operators still do not need Python or uv. Only the package preparation and returned-result comparison use uv.
