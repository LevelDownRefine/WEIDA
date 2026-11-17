"""Generate synthetic fixtures, run portable Octave, and report stage errors."""
import os
os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
import scipy
import sklearn
from scipy.io import savemat, loadmat
from scipy.signal import butter, detrend, filtfilt, hilbert
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from weida import compute_weighted_eigenvectors, analyze_state_occupancy, calculate_state_dynamics


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--octave',required=True)
    parser.add_argument('--output',default='validation_output')
    args=parser.parse_args()
    out=Path(args.output).resolve(); out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(284); tr=2.1
    time=np.arange(210)*tr
    offsets=rng.uniform(-np.pi,np.pi,(78,1))
    bold=np.sin(2*np.pi*.055*time+offsets)+.3*np.sin(2*np.pi*.025*time+offsets)+.15*rng.normal(size=(78,210))
    b,a=butter(2,[.04,.07],btype='bandpass',fs=1/tr)
    ts=detrend(bold,axis=1); ts-=ts.mean(axis=1,keepdims=True)
    filtered=filtfilt(b,a,ts,axis=1,padlen=3*(max(len(a),len(b))-1))
    phase=np.angle(hilbert(filtered,axis=1))
    # Four well-separated clusters: matching initial centers, not just RNG seeds.
    initial=np.eye(4)*8
    cluster_input=np.concatenate([center+rng.normal(scale=.1,size=(24,4)) for center in initial])
    sequences=[np.array([0,0,1,1,0,2,2,3,3,0]),np.array([2,2,2,2]),np.array([1,0,1,0,1,3,3]),np.array([3,3,0,0,2,1,1,0])]
    labels=np.concatenate(sequences); lengths=np.array(list(map(len,sequences))); conditions=np.array([1,1,2,2])
    savemat(out/'input.mat',dict(bold=bold,tr=tr,shared_phase=phase,cluster_input=cluster_input,initial_centers=initial,labels=(labels+1)[:,None],lengths=lengths[:,None],conditions=conditions[:,None]))
    reference=Path(__file__).with_name('reference.m')
    (out/'reference.m').write_bytes(reference.read_bytes())
    completed=subprocess.run([str(Path(args.octave).resolve()),'--no-gui','--quiet','reference.m'],cwd=out,capture_output=True,text=True)
    (out/'octave.log').write_text(completed.stdout+'\n'+completed.stderr,encoding='utf-8')
    if completed.returncode: raise RuntimeError(f'Octave failed; see {out / "octave.log"}')
    ref=loadmat(out/'reference.mat',squeeze_me=True)
    checks={}
    def check(name,left,right,tolerance=1e-9,required=True):
        left,right=np.asarray(left),np.asarray(right)
        error=float(np.max(np.abs(left[np.isfinite(left)&np.isfinite(right)]-right[np.isfinite(left)&np.isfinite(right)]),initial=0))
        checks[name]=dict(max_absolute_error=error,tolerance=tolerance,passed=bool(np.allclose(left,right,atol=tolerance,rtol=0,equal_nan=True)),required=required)
    check('filter_b',b,ref['b']); check('filter_a',a,ref['a'])
    check('filtered',filtered,ref['filtered'],1e-8,False)
    check('circular_phase',np.angle(np.exp(1j*(phase-ref['phase']))),np.zeros_like(phase),1e-8,False)
    features,_=compute_weighted_eigenvectors([bold],tr)
    check('canonical_features_shared_phase',features,ref['shared_features'],1e-8)
    check('canonical_features_end_to_end',features,ref['canonical_features'],1e-8,False)
    check('original_sign_features_end_to_end',features,ref['native_features'],1e-8,False)
    model=KMeans(n_clusters=4,init=initial,n_init=1,max_iter=500,algorithm='lloyd',random_state=0).fit(cluster_input)
    check('cluster_labels_fixed_start',model.labels_+1,ref['idx'],0)
    check('cluster_centers_fixed_start',model.cluster_centers_,ref['C'])
    distances=model.transform(cluster_input)**2
    check('cluster_distances',distances,ref['D'])
    check('cluster_sumd',np.bincount(model.labels_,weights=distances[np.arange(96),model.labels_]),ref['sumd'])
    check('silhouette',silhouette_samples(cluster_input,model.labels_,metric='sqeuclidean'),ref['ss'])
    p,lt=analyze_state_occupancy(labels,lengths,4,tr)
    check('occupancy',p,ref['P']); check('dwell_seconds',lt,ref['LT'])
    dynamics=calculate_state_dynamics(labels,np.eye(4),lengths,p,lt,conditions)
    check('state_order',dynamics['order']+1,ref['order'],0)
    check('scan_transitions',dynamics['transitions'],ref['transitions'])
    check('group_transitions',np.stack([dynamics['group_transitions'][str(i)] for i in [1,2]]),ref['group'])
    report=dict(runtime=str(ref['runtime_version']),numpy=np.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,fixture_seed=284,synthetic_only=True,matlab_equivalence=False,checks=checks)
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))
    if any(not c['passed'] and c['required'] for c in checks.values()): sys.exit(1)

if __name__=='__main__': main()
