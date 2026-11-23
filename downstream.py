"""下游分析：固定中心预测、状态离散度、功能分离度及配对统计。"""
import numpy as np
from scipy.stats import ttest_rel, false_discovery_control
from weida import compute_weighted_eigenvectors, analyze_state_occupancy, _sequences


def predict_states(vectors, centers, lengths, tr=2.1):
    """按最近中心预测，保持训练状态编号，不重新聚类或排序。

    无完整状态段的持续时间为 NaN；原 WEiDA_fix_cluster 返回 0。
    """
    x,c=np.asarray(vectors,float),np.asarray(centers,float)
    if x.ndim!=2 or c.ndim!=2 or not len(c) or x.shape[1]!=c.shape[1] or not np.isfinite(x).all() or not np.isfinite(c).all():
        raise ValueError('特征与中心必须有限且脑区维度一致')
    labels=np.argmin(np.sum((x[:,None,:]-c[None,:,:])**2,axis=2),axis=1)
    seqs=_sequences(labels,lengths,len(c))
    p,lt=analyze_state_occupancy(labels,lengths,len(c),tr)
    transition=np.zeros((len(seqs),len(c),len(c)))
    for i,seq in enumerate(seqs):
        np.add.at(transition[i],(seq[:-1],seq[1:]),1)
        totals=transition[i].sum(axis=1,keepdims=True)
        np.divide(transition[i],totals,out=transition[i],where=totals>0)
    return dict(states=labels,occupancy=p,lifetime=lt,transitions=transition)


def predict_bold(scans,centers,tr=2.1):
    """BOLD 与训练模型须采用相同脑区顺序、滤波及符号约定。"""
    vectors,lengths=compute_weighted_eigenvectors(scans,tr)
    return dict(vectors=vectors,lengths=lengths,**predict_states(vectors,centers,lengths,tr))


def state_metastability(vectors,labels,centers,lengths):
    """按原 MSanalysis 公式：sum((x-center)^2)/(脑区数*(帧数-1))。

    这是围绕固定中心的离散度，不是 Kuramoto 序参量标准差。
    状态少于两帧时返回 NaN。
    """
    x,c=np.asarray(vectors,float),np.asarray(centers,float)
    if x.ndim!=2 or c.ndim!=2 or x.shape[1]!=c.shape[1] or not x.shape[1] or not np.isfinite(x).all() or not np.isfinite(c).all():
        raise ValueError('特征与中心形状或数值无效')
    seqs=_sequences(labels,lengths,len(c))
    if len(x)!=len(labels): raise ValueError('标签与特征长度不一致')
    out=np.full((len(seqs),len(c)),np.nan); start=0
    for i,seq in enumerate(seqs):
        chunk=x[start:start+len(seq)]; start+=len(seq)
        for k,center in enumerate(c):
            selected=chunk[seq==k]; n=len(selected)
            if n>1: out[i,k]=np.sum((selected-center)**2)/(x.shape[1]*(n-1))
    return out


def functional_segregation(scan,center,clip=1e-7):
    """按中心正/非正分成两模块，对 Pearson r 做 Fisher z 后求内外差。

    输入为脑区×时间。完美相关裁剪到 ±(1-clip)，避免无穷；记录裁剪边数。
    模块均需至少两个脑区，节点指标排除自身边。
    """
    x,c=np.asarray(scan,float),np.asarray(center,float)
    if x.ndim!=2 or c.shape!=(x.shape[0],) or x.shape[1]<3 or not np.isfinite(x).all() or not np.isfinite(c).all() or np.any(np.ptp(x,axis=1)==0):
        raise ValueError('无效 BOLD 或中心')
    if not 0<clip<1: raise ValueError('裁剪参数需在 0 与 1 之间')
    positive=c>0
    if min(positive.sum(),(~positive).sum())<2: raise ValueError('每个模块至少需要两个脑区')
    corr=np.corrcoef(x); upper=np.triu(np.ones_like(corr,dtype=bool),1)
    clipped=int(np.sum(upper & (np.abs(corr)>1-clip)))
    z=np.arctanh(np.clip(corr,-1+clip,1-clip)); same=positive[:,None]==positive[None,:]
    within=float(z[upper & same].mean()); across=float(z[upper & ~same].mean())
    nodal=np.empty(len(c))
    for i in range(len(c)):
        inside=same[i].copy(); inside[i]=False
        nodal[i]=z[i,inside].mean()-z[i,~same[i]].mean()
    return dict(within=within,across=across,segregation=within-across,nodal=nodal,clipped_edges=clipped)


def paired_comparison(values,subject_ids,conditions,before='OFF',after='ON'):
    """按 ID 严格配对，逐特征成对排除 NaN，双侧 t 检验与 BH-FDR。

    每次调用的全部特征构成一个校正族；效应方向为 after-before。
    零差值方差时不作退化 t 推断，t/p/q 返回 NaN。
    """
    x=np.asarray(values,float); ids=list(subject_ids); groups=list(conditions)
    if x.ndim==1: x=x[:,None]
    if x.ndim!=2 or len(x)!=len(ids) or len(x)!=len(groups) or np.isinf(x).any() or before==after: raise ValueError('无效统计输入')
    maps=[]
    for group in [before,after]:
        mapping={}
        for i,(subject,condition) in enumerate(zip(ids,groups)):
            if condition==group:
                if subject in mapping: raise ValueError('重复受试者条件')
                mapping[subject]=i
        maps.append(mapping)
    if not maps[0] or maps[0].keys()!=maps[1].keys(): raise ValueError('配对不完整')
    paired=sorted(maps[0]); a=x[[maps[0][s] for s in paired]]; b=x[[maps[1][s] for s in paired]]
    result={key:np.full(x.shape[1],np.nan) for key in ['mean_before','mean_after','mean_difference','t','p','q']}
    result['n_pairs']=np.zeros(x.shape[1],dtype=int)
    for j in range(x.shape[1]):
        valid=np.isfinite(a[:,j]) & np.isfinite(b[:,j]); av=a[valid,j];bv=b[valid,j]
        n=len(av);result['n_pairs'][j]=n
        if n:
            result['mean_before'][j]=av.mean();result['mean_after'][j]=bv.mean();result['mean_difference'][j]=(bv-av).mean()
        if n>=2 and np.std(bv-av)>0:
            result['t'][j],result['p'][j]=ttest_rel(bv,av)
    finite=np.isfinite(result['p'])
    if finite.any(): result['q'][finite]=false_discovery_control(result['p'][finite],method='bh')
    result['paired_ids']=paired
    return result


def independent_comparison(values,conditions,reference='HC',target='OFF',equal_var=True):
    """独立样本双侧 t 检验；默认等方差以对应原脚本，可显式选择 Welch。"""
    from scipy.stats import ttest_ind
    x=np.asarray(values,float); groups=np.asarray(conditions)
    if x.ndim==1: x=x[:,None]
    if x.ndim!=2 or groups.shape!=(len(x),) or np.isinf(x).any() or reference==target: raise ValueError('无效统计输入')
    a=x[groups==reference]; b=x[groups==target]
    if not len(a) or not len(b): raise ValueError('缺少待比较组')
    result={key:np.full(x.shape[1],np.nan) for key in ['mean_reference','mean_target','mean_difference','t','p','q']}
    result['n_reference']=np.zeros(x.shape[1],int);result['n_target']=np.zeros(x.shape[1],int)
    for j in range(x.shape[1]):
        av=a[np.isfinite(a[:,j]),j];bv=b[np.isfinite(b[:,j]),j]
        result['n_reference'][j]=len(av);result['n_target'][j]=len(bv)
        if len(av) and len(bv):
            result['mean_reference'][j]=av.mean();result['mean_target'][j]=bv.mean();result['mean_difference'][j]=bv.mean()-av.mean()
        if len(av)>=2 and len(bv)>=2 and (np.std(av)>0 or np.std(bv)>0):
            result['t'][j],result['p'][j]=ttest_ind(bv,av,equal_var=equal_var)
    finite=np.isfinite(result['p'])
    if finite.any():result['q'][finite]=false_discovery_control(result['p'][finite],method='bh')
    return result
