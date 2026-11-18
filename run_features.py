"""从已有模型结果计算状态离散度、功能分离度及 STN 配对统计。"""
import argparse,csv,json
from pathlib import Path
import numpy as np
from downstream import state_metastability,functional_segregation,paired_comparison,predict_states,independent_comparison

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('input',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--module-state',type=int,default=2,help='用哪一个状态中心划分正/非正模块（从1开始）')
    args=p.parse_args();r=np.load(args.input/'results.npz');meta=json.loads((args.input/'metadata.json').read_text())
    k=len(r['centers']); state=args.module_state-1
    if not 0<=state<k: p.error('模块状态编号超出范围')
    args.output.mkdir(parents=True,exist_ok=False)
    ms=state_metastability(r['vectors'],r['states'],r['centers'],r['lengths'])
    fs=[];nodal=[];clip=[]
    for m in meta['scans']:
        with Path(m['file']).open(encoding='utf-8-sig',newline='') as f: rows=list(csv.reader(f))
        if rows[0][1:]!=meta['regions']: raise ValueError('脑区顺序与训练模型不一致')
        x=np.array([[float(v) for v in row[1:]] for row in rows[1:]]).T
        v=functional_segregation(x,r['centers'][state]); fs.append([v['within'],v['across'],v['segregation']]);nodal.append(v['nodal']);clip.append(v['clipped_edges'])
    fs=np.asarray(fs);nodal=np.asarray(nodal)
    ids=[m['subject'].removeprefix('sub-STNon').removeprefix('sub-STNoff') for m in meta['scans']]
    groups=[m['condition'] for m in meta['scans']]
    features={'occupancy':r['occupancy'],'lifetime':r['lifetime'],'metastability':ms,'functional_segregation':fs[:,2],'nodal_segregation':nodal}
    # 没有有效源状态的转移行不作为 0 纳入统计。
    trans=np.where(r['transition_valid'][:,:,None],r['transitions'],np.nan)
    features['transitions']=trans.reshape(len(ids),-1)
    stats={}
    for name,values in features.items():
        result=paired_comparison(values,ids,groups)
        stats[name]={key:(value.tolist() if isinstance(value,np.ndarray) else value) for key,value in result.items()}
    independent={}
    for target in ['OFF','ON']:
        independent[target+'_vs_HC']={}
        for name,values in features.items():
            result=independent_comparison(values,groups,target=target)
            independent[target+'_vs_HC'][name]={key:value.tolist() for key,value in result.items()}
    (args.output/'independent_statistics.json').write_text(json.dumps(independent,indent=2),encoding='utf-8')
    prediction=predict_states(r['vectors'],r['centers'],r['lengths'],meta['tr'])
    agreement=float(np.mean(prediction['states']==r['states']))
    np.savez_compressed(args.output/'features.npz',metastability=ms,functional_segregation=fs,nodal_segregation=nodal,clipped_edges=clip)
    (args.output/'paired_statistics.json').write_text(json.dumps(stats,indent=2),encoding='utf-8')
    info={'scans':meta['scans'],'module_state':args.module_state,'module_sizes':[int(np.sum(r['centers'][state]>0)),int(np.sum(r['centers'][state]<=0))],'fixed_center_training_agreement':agreement,'correction':'BH-FDR separately per feature family, ON minus OFF','notes':'探索性；同一数据训练后统计；模块状态编号不保证对应论文状态。无MATLAB等价验证。'}
    (args.output/'metadata.json').write_text(json.dumps(info,indent=2,ensure_ascii=False),encoding='utf-8')
    print('扫描数:',len(ids),'训练集固定中心标签一致率:',agreement,'结果:',args.output)

if __name__=='__main__':main()
