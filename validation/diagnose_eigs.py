"""Inspect raw solver behavior on byte-identical MATLAB matrices; never align signs."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.sparse.linalg import eigsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from weida import weight_eigenvectors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder')
    args = parser.parse_args()
    folder = Path(args.folder)
    ref = loadmat(folder / 'eigs_diagnosis.mat', simplify_cells=True)
    matrices = ref['matrices']
    target = np.array([weight_eigenvectors(w,v) for w,v in zip(ref['native_values'],ref['native_vectors'])])
    results = {}

    def compare(name, actual):
        error = np.max(np.abs(np.asarray(actual)-target),axis=1)
        results[name] = {'matching_cases': int((error < 1e-8).sum()),
                         'total_cases': len(target), 'max_absolute_error': float(error.max())}

    for mode in ('repeat','seed','handle','start'):
        compare('matlab_'+mode,[weight_eigenvectors(w,v)
                               for w,v in zip(ref[mode+'_values'],ref[mode+'_vectors'])])
    dense, iterative, residuals, correlations, eigenvalue_errors = [], [], [], [], []
    for i, matrix in enumerate(matrices):
        w,v = np.linalg.eigh(matrix)
        idx = np.argsort(np.abs(w))[-6:][::-1]
        w,v = w[idx],v[:,idx]
        dense.append(weight_eigenvectors(w,v))
        residuals.append(float(np.max(np.abs(matrix@v-v*w))))
        eigenvalue_errors.append(float(np.max(np.abs(w-ref['native_values'][i]))))
        correlations.append(np.diag(v[:,:2].T@ref['native_vectors'][i,:,:2]).tolist())
        trace = ref['operator_inputs'][i]
        w,v = eigsh(matrix,k=6,which='LM',v0=trace[:,0],ncv=20,tol=1e-14,maxiter=300)
        idx = np.argsort(np.abs(w))[::-1]
        iterative.append(weight_eigenvectors(w[idx],v[:,idx]))
    compare('numpy_eigh_identical_matrix',dense)
    compare('scipy_eigsh_matlab_first_operator_input',iterative)
    report = {'runtime': str(ref['runtime_version']), 'all_converged': bool((ref['flags']==0).all()),
              'results': results, 'numpy_max_eigenpair_residual': max(residuals),
              'numpy_max_eigenvalue_error': max(eigenvalue_errors),
              'dominant_signed_vector_correlations': correlations,
              'operator_call_counts': [int(x.shape[1]) for x in ref['operator_inputs']],
              'profile_names': np.atleast_1d(ref['profile_names']).tolist()}
    (folder / 'eigs_diagnosis.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('profile_names','dominant_signed_vector_correlations')},indent=2))


if __name__ == '__main__':
    main()
