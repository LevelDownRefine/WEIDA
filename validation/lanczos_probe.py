"""Diagnostic Lanczos probe, not a complete MATLAB eigs replacement."""
import numpy as np
from scipy.linalg import eigh_tridiagonal


def lanczos_probe(matrix, start, steps=20):
    """Stop on loss of numerical independence; never divide by that residual.

    The relative threshold uses machine precision, not the reference features.
    This stops the experiment rather than implementing Krylov-Schur restart.
    """
    matrix = np.asarray(matrix, dtype=float)
    steps = min(steps, len(matrix))
    basis = np.zeros((len(matrix), steps))
    basis[:, 0] = start
    diagonal = np.zeros(steps)
    off_diagonal = np.zeros(steps - 1)
    breakdown = False
    for j in range(steps):
        residual = matrix @ basis[:, j]
        diagonal[j] = basis[:, j] @ residual
        residual -= diagonal[j] * basis[:, j]
        if j:
            residual -= off_diagonal[j - 1] * basis[:, j - 1]
        before = np.linalg.norm(residual)
        for _ in range(2):
            residual -= basis[:, :j + 1] @ (basis[:, :j + 1].T @ residual)
        beta = np.linalg.norm(residual)
        breakdown = beta <= np.finfo(float).eps * before
        if breakdown or j == steps - 1:
            break
        off_diagonal[j] = beta
        basis[:, j + 1] = residual / beta
    values, vectors = eigh_tridiagonal(
        diagonal[:j + 1], off_diagonal[:j], lapack_driver='stev')
    index = np.argsort(np.abs(values))[-6:][::-1]
    return values[index], basis[:, :j + 1] @ vectors[:, index], {
        'steps': j + 1, 'breakdown': bool(breakdown)}
