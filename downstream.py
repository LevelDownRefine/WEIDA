"""下游分析：固定中心预测、状态离散度、功能分离度及配对统计。"""
import numpy as np
from research_core.statistics import paired_comparison, independent_comparison
from research_core.networks import module_summary
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
    z=np.arctanh(np.clip(corr,-1+clip,1-clip))
    return dict(**module_summary(z, positive, directed=False), clipped_edges=clipped)
