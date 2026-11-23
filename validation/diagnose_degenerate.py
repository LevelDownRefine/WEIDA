"""Locate numerical amplification and distinguish basis from subspace error."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from validation.lanczos_probe import lanczos_probe
from validation.matlab_start_probe import default_start
from weida import weight_eigenvectors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder')
    args = parser.parse_args()
    folder = Path(args.folder)
    ref = loadmat(folder / 'eigs_diagnosis.mat', simplify_cells=True)
    # Explicitly fix memory layout so local traces can be reproduced.
    matrix = np.array(ref['matrices'][-1], order='C')
    expected_vectors = ref['native_vectors'][-1]
    expected_values = ref['native_values'][-1]
    expected = weight_eigenvectors(expected_values, expected_vectors)
    reference_basis = ref['operator_inputs'][-1][:, 1:]
    results = {}
    for name, start in [('shared', reference_basis[:, 0]),
                        ('independent', default_start(len(matrix)))]:
        basis = np.zeros((len(matrix), 6))
        basis[:, 0] = start
        beta = 0.
        trace = []
        for j in range(5):
            residual = matrix @ basis[:, j]
            alpha = basis[:, j] @ residual
            residual -= alpha * basis[:, j]
            if j:
                residual -= beta * basis[:, j - 1]
            for _ in range(2):
                residual -= basis[:, :j + 1] @ (basis[:, :j + 1].T @ residual)
            beta = np.linalg.norm(residual)
            if beta == 0:
                break
            basis[:, j + 1] = residual / beta
            trace.append({'iteration': j + 1, 'residual_norm': float(beta),
                          'next_basis_max_error': float(np.max(np.abs(
                              basis[:, j + 1] - reference_basis[:, j + 1])))})
        values, vectors, _ = lanczos_probe(matrix, start)
        actual = weight_eigenvectors(values, vectors)
        projector_error = vectors[:, :2] @ vectors[:, :2].T
        projector_error -= expected_vectors[:, :2] @ expected_vectors[:, :2].T
        results[name] = {
            'trace': trace,
            'feature_max_error': float(np.max(np.abs(actual - expected))),
            'feature_error_allowing_global_flip': float(min(
                np.max(np.abs(actual - expected)), np.max(np.abs(actual + expected)))),
            'projector_max_error': float(np.max(np.abs(projector_error))),
            'eigenpair_residual': float(np.max(np.abs(matrix @ vectors - vectors * values))),
            'dominant_basis_overlap': (vectors[:, :2].T @ expected_vectors[:, :2]).tolist()}
    report = {'top_eigenvalues': expected_values[:2].tolist(),
              'top_gap': float(abs(expected_values[0] - expected_values[1])),
              'numpy_version': np.__version__, 'matrix_layout': 'C', 'results': results}
    (folder / 'degenerate_python.json').write_text(json.dumps(report, indent=2),encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
