"""Generate shared fixtures and compare MATLAB/Octave with production Python."""
import os

os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
import sklearn
from scipy.io import loadmat, savemat
from scipy.signal import butter, hilbert
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from preprocessing import preprocess_bold
from weida import (
    analyze_state_occupancy,
    calculate_state_dynamics,
    compute_weighted_eigenvectors,
    weight_eigenvectors,
)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--octave')
    parser.add_argument('--prepare-only',action='store_true',help='Export a standalone MATLAB/Octave handoff')
    parser.add_argument('--compare-only',action='store_true',help='Compare returned reference.mat without an M runtime')
    parser.add_argument('--output',default='validation_output')
    parser.add_argument('--strict-signal',action='store_true',help='Require filtered signals and phases to match (original features are always required)')
    args=parser.parse_args()
    if sum([bool(args.octave),args.prepare_only,args.compare_only]) != 1:
        parser.error('Choose exactly one of --octave, --prepare-only, --compare-only')
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(284); tr=2.1
    time=np.arange(210)*tr
    offsets=rng.uniform(-np.pi,np.pi,(78,1))
    bold=np.sin(2*np.pi*.055*time+offsets)+.3*np.sin(2*np.pi*.025*time+offsets)+.15*rng.normal(size=(78,210))
    b,a=butter(2,[.04,.07],btype='bandpass',fs=1/tr)
    filtered=preprocess_bold(bold,tr=tr)
    phase=np.angle(hilbert(filtered,axis=1))
    # Four well-separated clusters: matching initial centers, not just RNG seeds.
    initial=np.eye(4)*8
    cluster_input=np.concatenate([center+rng.normal(scale=.1,size=(24,4)) for center in initial])
    sequences=[np.array([0,0,1,1,0,2,2,3,3,0]),np.array([2,2,2,2]),np.array([1,0,1,0,1,3,3]),np.array([3,3,0,0,2,1,1,0])]
    labels=np.concatenate(sequences); lengths=np.array(list(map(len,sequences))); conditions=np.array([1,1,2,2])
    if not args.compare_only:
        savemat(out/'input.mat',{'bold': bold,'tr': tr,'shared_phase': phase,'shared_filtered': filtered,'cluster_input': cluster_input,'initial_centers': initial,'labels': (labels+1)[:,None],'lengths': lengths[:,None],'conditions': conditions[:,None]})
        for name in ['reference.m','run_reference.m','HANDOFF.md']:
            (out/name).write_bytes(Path(__file__).with_name(name).read_bytes())
        manifest={name:hashlib.sha256((out/name).read_bytes()).hexdigest() for name in ['input.mat','reference.m','run_reference.m']}
        (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    manifest=json.loads((out/'manifest.json').read_text(encoding='utf-8'))
    for name,digest in manifest.items():
        if hashlib.sha256((out/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(f'Handoff checksum mismatch: {name}')
    if args.prepare_only:
        print(f'Handoff ready: {out}'); return
    if args.octave:
        completed=subprocess.run([str(Path(args.octave).resolve()),'--no-gui','--quiet','--eval','run_reference'],cwd=out,capture_output=True,text=True,check=False)
        (out/'octave.log').write_text(completed.stdout+'\n'+completed.stderr,encoding='utf-8')
        if completed.returncode:
            print(completed.stdout, file=sys.stderr); print(completed.stderr, file=sys.stderr)
            raise RuntimeError(f'Octave failed; see {out / "octave.log"}')
    # Always compare against the actual handed-off fixture, not a regenerated RNG stream.
    fixture=loadmat(out/'input.mat',squeeze_me=True)
    bold=fixture['bold']; tr=float(fixture['tr']); phase=fixture['shared_phase']
    cluster_input=fixture['cluster_input']; initial=fixture['initial_centers']
    labels=fixture['labels'].astype(int)-1; lengths=fixture['lengths'].astype(int); conditions=fixture['conditions'].astype(int)
    b,a=butter(2,[.04,.07],btype='bandpass',fs=1/tr)
    filtered=preprocess_bold(bold,tr=tr)
    # Older handoff packages did not include shared_filtered.
    shared_features,_=compute_weighted_eigenvectors([fixture.get('shared_filtered',filtered)],tr,preprocessed=True)
    phase=np.angle(hilbert(filtered,axis=1))
    ref=loadmat(out/'reference.mat',squeeze_me=True)
    checks={}
    def check(name,left,right,tolerance=1e-9,required=True):
        left,right=np.asarray(left),np.asarray(right)
        if left.shape != right.shape:
            checks[name]={'passed': False,'required': required,'error': f'Shape mismatch: {left.shape} vs {right.shape}'}; return
        error=float(np.max(np.abs(left[np.isfinite(left)&np.isfinite(right)]-right[np.isfinite(left)&np.isfinite(right)]),initial=0))
        checks[name]={'max_absolute_error': error,'tolerance': tolerance,'passed': bool(np.allclose(left,right,atol=tolerance,rtol=0,equal_nan=True)),'required': required}
    check('filter_b',b,ref['b']); check('filter_a',a,ref['a'])
    check('filtered',filtered,ref['filtered'],1e-8,args.strict_signal)
    check('circular_phase',np.angle(np.exp(1j*(phase-ref['phase']))),np.zeros_like(phase),1e-8,args.strict_signal)
    weighted=np.array([weight_eigenvectors(values,vectors)
                       for values,vectors in zip(ref['native_eigenvalues'],ref['native_eigenvectors'])])
    check('weighted_features_shared_eigenpairs',weighted,ref['native_features'],1e-8)
    features,_=compute_weighted_eigenvectors([bold],tr)
    check('original_features_shared_phase',shared_features,ref['shared_features'],1e-8)
    check('original_features_end_to_end',features,ref['native_features'],1e-8)
    model=KMeans(n_clusters=4,init=initial,n_init=1,max_iter=500,algorithm='lloyd',random_state=0).fit(cluster_input)
    check('cluster_labels_fixed_start',model.labels_+1,ref['idx'],0)
    check('cluster_centers_fixed_start',model.cluster_centers_,ref['C'])
    distances=model.transform(cluster_input)**2
    check('cluster_distances',distances,ref['D'])
    check('cluster_sumd',np.bincount(model.labels_,weights=distances[np.arange(len(cluster_input)),model.labels_]),ref['sumd'])
    check('silhouette',silhouette_samples(cluster_input,model.labels_,metric='sqeuclidean'),ref['ss'])
    p,lt=analyze_state_occupancy(labels,lengths,4,tr)
    check('occupancy',p,ref['P']); check('dwell_seconds',lt,ref['LT'])
    dynamics=calculate_state_dynamics(labels,np.eye(4),lengths,p,lt,conditions)
    check('state_order',dynamics['order']+1,ref['order'],0)
    check('scan_transitions',dynamics['transitions'],ref['transitions'])
    check('group_transitions',np.stack([dynamics['group_transitions'][str(i)] for i in [1,2]]),ref['group'])
    report={'input_sha256': manifest['input.mat'],'runtime': str(ref['runtime_version']),'numpy': np.__version__,'scipy': scipy.__version__,'sklearn': sklearn.__version__,'fixture_seed': 284,'synthetic_only': True,'matlab_equivalence': False,'strict_signal': args.strict_signal,'checks': checks}
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    if any(not c['passed'] and c['required'] for c in checks.values()): sys.exit(1)

if __name__=='__main__': main()
