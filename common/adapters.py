"""将 WEIDA 已有结果转换为共享指标表，不重新拟合状态。"""
import json
from pathlib import Path

import numpy as np

from .postprocessing import FeatureTable


def load_weida(directory,features_directory=None):
    directory=Path(directory)
    meta=json.loads((directory/'metadata.json').read_text(encoding='utf-8'))
    arrays=[];names=[];families=[];units=[]
    def add(family,values,unit,labels):
        values=np.asarray(values,float).reshape(len(meta['scans']),-1)
        if values.shape[1]!=len(labels): raise ValueError('WEiDA feature dimension mismatch')
        arrays.append(values);names.extend(labels);families.extend([family]*len(labels));units.extend([unit]*len(labels))
    with np.load(directory/'results.npz') as z:
        k=z['occupancy'].shape[1]
        if z['lifetime'].shape!=z['occupancy'].shape or z['transitions'].shape!=(len(meta['scans']),k,k):
            raise ValueError('WEiDA state axes do not match')
        add('occupancy',z['occupancy'],'fraction',[f'FO_state{i+1}' for i in range(k)])
        add('lifetime',z['lifetime'],'s',[f'DT_state{i+1}' for i in range(k)])
        trans=np.where(z['transition_valid'][:,:,None],z['transitions'],np.nan)
        add('transitions',trans,'probability',[f'TP_{i+1}_to_{j+1}' for i in range(k) for j in range(k)])
    if features_directory is not None:
        fd=Path(features_directory)
        feature_meta=json.loads((fd/'metadata.json').read_text(encoding='utf-8'))
        if [(m['subject'],m['condition']) for m in feature_meta['scans']]!=[(m['subject'],m['condition']) for m in meta['scans']]:
            raise ValueError('WEiDA downstream scan order differs from model result')
        with np.load(fd/'features.npz') as z:
            add('metastability',z['metastability'],'WEiDA variance',[f'VAR_state{i+1}' for i in range(k)])
            # 上游输出列为模块内、模块间和差值；这里只选择分离度差值。
            add('functional_segregation',z['functional_segregation'][:,2],'Fisher-z difference',['FS'])
            add('nodal_segregation',z['nodal_segregation'],'Fisher-z difference',[f'FS_{r}' for r in meta['regions']])
    return FeatureTable(np.concatenate(arrays,axis=1),meta['scans'],tuple(names),tuple(families),tuple(units),
        'WEiDA','Existing WEIDA Python model outputs; fixed state labels; no refitting')
