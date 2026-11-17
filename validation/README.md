# Octave comparison harness

`compare.py` creates a deterministic synthetic 78-region, 210-frame dataset, runs `reference.m` in Octave, and writes a structured `report.json` and intermediate `.mat` files. No clinical data is uploaded or required.

```
uv run --locked python validation/compare.py --octave C:/path/to/octave-cli.exe --output validation_output
```

The reference implements the source equations independently in M syntax, with Octave signal/statistics packages. It is not the complete untouched upstream script: data loading and interactive plotting are replaced with deterministic fixtures; source eigenvector-sign behavior and canonicalized behavior are both exported. K-means is tested with shared initial centers on separated clusters (not the source's random 300 restarts). State metrics use prescribed labels, including missing states, constant-state scans and incomplete boundary runs.

Checks distinguish:

- Filter coefficients and Octave vs SciPy filtering/phase differences.
- Original-sign and canonicalized weighted features, including a shared-phase comparison isolating the eigensolver.
- Fixed-initialization K-means centers, labels, squared distances, SUMD and silhouette.
- Occupancy, dwell time, sorted state order, individual transition rows and group means.

Filtering and original-sign end-to-end differences are reported without forcing equality. The executable exits nonzero when a required isolated-stage check fails. Passing does not establish MATLAB equivalence or scientific validity. Numeric reports are tied to the recorded Octave version and generated fixtures.

## Recorded local verification

Portable Octave 11.3.0 with signal 1.4.7 and statistics 1.8.3 was tested on Windows. See `report_octave_11.3.0.json` for exact errors and Python library versions. Of 17 checks, 16 pass; original-sign features intentionally remain a documented mismatch. The other checks include filtering/phase, canonicalized eigenfeatures, shared-initialization clustering and prescribed-label state statistics. The mismatch is evidence that the current canonicalization is a method change, not a faithful reproduction of the source's eigensolver sign choices. This test does not run the full original script or validate real-data results.

## Unit tests, CI and MATLAB handoff

Run `uv run --locked python -m unittest discover -v` from the repository root. Unit tests use hand-calculated expectations (not snapshots of the implementation): boundary/absent-state dwell, scan boundary isolation, scan-weighted group averages, state remapping and invalid input rejection. A handoff test verifies package generation without an M runtime and checks detection of altered scripts.

GitHub Actions runs Python tests and a separate Octave comparison on pushes and PRs. The Octave job installs signal/statistics from Ubuntu packages and uploads intermediate results, logs and report even on failure. Its green status only asserts the required isolated-stage checks; the original-sign feature difference remains reported, not declared equivalent. Runtime versions are recorded because Ubuntu Octave may differ from the local portable runtime.

Prepare a handoff without MATLAB/Octave:

```
uv run --locked python validation/compare.py --prepare-only --output validation_output/matlab_handoff
```

Zip that directory and give it to the MATLAB operator. `HANDOFF.md` contains Chinese instructions. They run `run_reference` and return the folder with `reference.mat` and `runtime.log`. Input/script hashes and version logs accompany the result. The same handoff can be rehearsed with Octave.

Validate returned results without MATLAB/Octave:

```
uv run --locked python validation/compare.py --compare-only --output /path/to/returned-folder
```

Retain the original handoff manifest locally and compare it to the returned copy; the bundled manifest detects accidental edits but is not a signature against deliberate replacement. Checks cover the independent reference harness, not all original upstream code or actual clinical data. MATLAB execution remains pending until an operator returns results.
