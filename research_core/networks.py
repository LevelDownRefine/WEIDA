"""Module summaries for matrices with explicitly declared edge semantics."""
import numpy as np


def module_summary(matrix, modules, *, directed):
    """Summarize within/across-module weights, excluding all self edges.

    Undirected matrices must be symmetric and use each upper-triangle edge
    once. Directed matrices use BOTH orientations; nodal values summarize
    outgoing edges. No Fisher transform or implicit X/Y symmetrization.
    NaN on the diagonal is allowed, but missing off-diagonal values are not.
    """
    x = np.asarray(matrix, dtype=float)
    labels = np.asarray(modules)
    if x.ndim != 2 or x.shape[0] != x.shape[1] or labels.shape != (len(x),):
        raise ValueError('Expected square matrix and one module label per region')
    off = ~np.eye(len(x), dtype=bool)
    if not np.isfinite(x[off]).all():
        raise ValueError('Nonfinite off-diagonal weights')
    same = labels[:, None] == labels[None, :]
    inside = same & off
    outside = ~same & off
    if not inside.any(axis=1).all() or not outside.any(axis=1).all():
        raise ValueError('Each region needs within- and across-module edges')
    if not directed and not np.allclose(x[off], x.T[off], rtol=1e-12, atol=1e-14):
        raise ValueError('Undirected matrix must be symmetric')
    edges = off if directed else np.triu(off, 1)
    within = float(x[edges & same].mean())
    across = float(x[edges & ~same].mean())
    nodal = np.array([x[i, inside[i]].mean() - x[i, outside[i]].mean()
                      for i in range(len(x))])
    return dict(within=within, across=across, segregation=within-across, nodal=nodal)
