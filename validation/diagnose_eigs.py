"""Inspect raw solver behavior on byte-identical MATLAB matrices; never align signs."""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.signal import hilbert
from scipy.sparse.linalg import eigsh

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from preprocessing import preprocess_bold
from validation.lanczos_probe import lanczos_probe
from validation.matlab_start_probe import DefaultStream, default_start
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
        if 'case_groups' in ref:
            groups = np.asarray(ref['case_groups'])
            results[name]['groups'] = {
                group: {'matching_cases': int((error[groups == group] < 1e-8).sum()),
                        'total_cases': int((groups == group).sum()),
                        'max_absolute_error': float(error[groups == group].max())}
                for group in dict.fromkeys(groups.tolist())}

    for mode in ('repeat','seed','handle','start'):
        compare('matlab_'+mode,[weight_eigenvectors(w,v)
                               for w,v in zip(ref[mode+'_values'],ref[mode+'_vectors'])])
    dense, iterative, lanczos, residuals, correlations, eigenvalue_errors = [], [], [], [], [], []
    lanczos_details = []
    independent = []
    generated_start = default_start(len(matrices[0]))
    for i, matrix in enumerate(matrices):
        w,v = np.linalg.eigh(matrix)
        idx = np.argsort(np.abs(w))[-6:][::-1]
        w,v = w[idx],v[:,idx]
        dense.append(weight_eigenvectors(w,v))
        residuals.append(float(np.max(np.abs(matrix@v-v*w))))
        eigenvalue_errors.append(float(np.max(np.abs(w-ref['native_values'][i]))))
        correlations.append(np.diag(v[:,:2].T@ref['native_vectors'][i,:,:2]).tolist())
        trace = ref['operator_inputs'][i]
        w,v = eigsh(matrix,k=6,which='LM',v0=trace[:,1],ncv=20,tol=1e-14,maxiter=300)
        idx = np.argsort(np.abs(w))[::-1]
        iterative.append(weight_eigenvectors(w[idx],v[:,idx]))
        # Diagnostic reconstruction: MATLAB's starting vector, then Python-only
        # Lanczos, two-pass reorthogonalization and tridiagonal QR. No sign matching.
        w,v,details = lanczos_probe(matrix,trace[:,1])
        lanczos.append(weight_eigenvectors(w,v))
        lanczos_details.append(details)
        w,v,_ = lanczos_probe(matrix,generated_start)
        independent.append(weight_eigenvectors(w,v))
    compare('numpy_eigh_identical_matrix',dense)
    compare('scipy_eigsh_matlab_start',iterative)
    compare('python_lanczos_matlab_start',lanczos)
    compare('python_lanczos_independent_start',independent)
    rng_observations = None
    if 'rng_fixture' in ref:
        fixture = ref['rng_fixture']
        stream = DefaultStream()
        uniforms = np.array([stream.uniform() for _ in fixture['uniform']])
        normals = DefaultStream().normal(len(fixture['normal']))
        rng_observations = {
            'sample_count': len(normals),
            'uniform_max_error': float(np.max(np.abs(uniforms-fixture['uniform']))),
            'normal_max_error': float(np.max(np.abs(normals-fixture['normal']))),
            'start_max_error': float(np.max(np.abs(generated_start-ref['operator_inputs'][0][:,1])))}
        np.testing.assert_array_equal(uniforms,fixture['uniform'])
        np.testing.assert_allclose(normals,fixture['normal'],rtol=0,atol=1e-12)
    octave_path = folder / 'octave_eigs_diagnosis.mat'
    end_to_end = None
    if (folder / 'input.mat').exists() and (folder / 'reference.mat').exists():
        inputs = loadmat(folder / 'input.mat', simplify_cells=True)
        expected = loadmat(folder / 'reference.mat', simplify_cells=True)['native_features']
        phase = np.angle(hilbert(preprocess_bold(inputs['bold'],tr=float(inputs['tr'])),axis=1))
        actual = []
        for phi in phase.T:
            w,v,_ = lanczos_probe(np.cos(phi[:,None]-phi[None,:]),generated_start)
            actual.append(weight_eigenvectors(w,v))
        errors = np.max(np.abs(np.array(actual)-expected),axis=1)
        end_to_end = {'matching_frames': int((errors < 1e-8).sum()),
                      'total_frames': len(errors), 'max_absolute_error': float(errors.max())}
    octave_summary = None
    if octave_path.exists():
        octave = loadmat(octave_path, simplify_cells=True)
        compare('octave_identical_matlab_matrix', octave['features'])
        octave_summary = {'runtime': str(octave['runtime_version']),
                          'all_converged': bool((octave['flags'] == 0).all())}
    report = {'runtime': str(ref['runtime_version']), 'all_converged': bool((ref['flags']==0).all()),
              'octave': octave_summary,
              'lanczos_details': lanczos_details,
              'independent_rng': rng_observations,
              'independent_end_to_end': end_to_end,
              'results': results, 'numpy_max_eigenpair_residual': max(residuals),
              'numpy_max_eigenvalue_error': max(eigenvalue_errors),
              'dominant_signed_vector_correlations': correlations,
              'operator_call_counts': [int(x.shape[1]) for x in ref['operator_inputs']],
              'profile_names': np.atleast_1d(ref['profile_names']).tolist(),
              'closest_stream_trials': sorted(ref['stream_trials'],key=lambda x:x['error'])[:5]}
    (folder / 'eigs_diagnosis.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('profile_names','dominant_signed_vector_correlations','lanczos_details')},indent=2))


if __name__ == '__main__':
    main()
