"""Two-sided study comparisons; independent of the feature-producing algorithm."""
import numpy as np
from scipy.stats import ttest_rel, false_discovery_control


def paired_comparison(values,subject_ids,conditions,before='OFF',after='ON'):
    """按 ID 严格配对，逐特征成对排除 NaN，双侧 t 检验与 BH-FDR。

    每次调用中有限的 p 值构成一个校正族；效应方向为 after-before。
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


def compare_feature_families(features, subject_ids, conditions):
    """Common study entry: {family: scan × feature array} -> three contrasts.

    Rows follow subject_ids/conditions; patient IDs must match across ON/OFF.
    Each family and contrast has its own BH correction over finite p values,
    preserving the existing WEiDA policy. Never pool edges as participants.
    """
    return {
        'ON_vs_OFF': {name: paired_comparison(values, subject_ids, conditions)
                      for name, values in features.items()},
        **{target + '_vs_HC': {
            name: independent_comparison(values, conditions, target=target)
            for name, values in features.items()}
            for target in ('OFF', 'ON')},
    }


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
