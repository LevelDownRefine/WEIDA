"""与算法无关的受试者指标表、统计检验、作图及结果保存。"""
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy import stats


@dataclass
class FeatureTable:
    values: np.ndarray
    scans: list
    names: tuple
    families: tuple
    units: tuple
    method: str
    definition: str

    def __post_init__(self):
        self.values=np.asarray(self.values,float)
        n=len(self.names)
        if self.values.shape!=(len(self.scans),n) or n==0 or np.isinf(self.values).any():
            raise ValueError('Expected scans x features, permitting NaN but not infinity')
        if len(set(self.names))!=n or len(self.families)!=n or len(self.units)!=n:
            raise ValueError('Feature names must be unique; every feature needs family and unit')
        if not self.method or not self.definition or not all(self.families):
            raise ValueError('Method, feature definition and correction families are required')
        group_indices(self.scans)


def group_indices(scans):
    groups={g:{} for g in ('HC','OFF','ON')}
    prefixes={'HC':'sub-HC','OFF':'sub-STNoff','ON':'sub-STNon'}
    for i,m in enumerate(scans):
        g=m['condition']
        if g not in groups: raise ValueError(f'Unsupported condition: {g}')
        if 'participant_id' in m:
            participant=str(m['participant_id'])
        else:
            if not m['subject'].startswith(prefixes[g]):
                raise ValueError('Supply participant_id when subject does not use DBS naming')
            participant=m['subject'][len(prefixes[g]):]
        if not participant or participant in groups[g]:
            raise ValueError('Empty or duplicate participant-condition')
        groups[g][participant]=i
    if any(len(v)<2 for v in groups.values()):
        raise ValueError('Require at least two participants per condition')
    if groups['OFF'].keys()!=groups['ON'].keys():
        raise ValueError('Incomplete ON/OFF pairs')
    return {g:[v[k] for k in sorted(v)] for g,v in groups.items()}


def adjust_bh(p):
    p=np.asarray(p,float);valid=np.isfinite(p)
    q=stats.false_discovery_control(np.where(valid,p,1.).ravel()).reshape(p.shape)
    return np.where(valid,q,np.nan)


def compare_features(table):
    """双侧 t 检验；每个对比 × 声明的指标族分别进行 BH 校正。

    配对缺失值成对排除，独立组按指标排除。
    不可估计项保留预定校正族中的位置，校正时按 p=1 处理。
    """
    idx=group_indices(table.scans);rows=[]
    for first,second in [('ON','OFF'),('OFF','HC'),('ON','HC')]:
        current=[];paired=second=='OFF'
        for k,name in enumerate(table.names):
            a,b=table.values[idx[first],k],table.values[idx[second],k]
            if paired:
                good=np.isfinite(a)&np.isfinite(b);a,b=a[good],b[good]
            else:
                a,b=a[np.isfinite(a)],b[np.isfinite(b)]
            difference=float(a.mean()-b.mean()) if len(a) and len(b) else np.nan
            t=p=low=high=np.nan
            if min(len(a),len(b))>=2:
                df=len(a)-1 if paired else len(a)+len(b)-2
                se=((a-b).std(ddof=1)/np.sqrt(len(a)) if paired else
                    np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/df*(1/len(a)+1/len(b))))
                if se>0 and np.isfinite(se):
                    test=stats.ttest_rel(a,b) if paired else stats.ttest_ind(a,b,equal_var=True)
                    t,p=float(test.statistic),float(test.pvalue)
                    margin=stats.t.ppf(.975,df)*se;low,high=difference-margin,difference+margin
            current.append({'contrast': first+'-'+second,'feature': name,'family': table.families[k],
                'unit': table.units[k],'n_first': len(a),'n_second': len(b),'difference': difference,
                'ci95_low': low,'ci95_high': high,'t': t,'p': p,'q': np.nan,
                'test': 'paired_t' if paired else 'independent_equal_variance_t'})
        for family in dict.fromkeys(table.families):
            selected=[r for r in current if r['family']==family]
            for r,q in zip(selected,adjust_bh([r['p'] for r in selected])): r['q']=float(q)
        rows.extend(current)
    return rows


def plot_features(table,output):
    """为两种算法使用相同的分布图和患者配对变化图。"""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    idx=group_indices(table.scans);colors=['#7c899b','#d58b40','#5b9d9a']
    # 按校正族分别作图，避免将不同量纲的指标混在一起。
    families=list(dict.fromkeys(table.families))
    for family in families:
        selected=[k for k,f in enumerate(table.families) if f==family]
        for start in range(0,len(selected),16):
            subset=selected[start:start+16];ncols=min(4,len(subset));nrows=(len(subset)+ncols-1)//ncols
            for paired,kind in [(False,'group_distributions'),(True,'paired_changes')]:
                fig,axes=plt.subplots(nrows,ncols,figsize=(3.3*ncols,2.8*nrows),squeeze=False)
                for k,ax in zip(subset,axes.flat):
                    if paired:
                        for a,b in zip(table.values[idx['OFF'],k],table.values[idx['ON'],k]):
                            if np.isfinite([a,b]).all(): ax.plot([0,1],[a,b],color=colors[2],alpha=.5,lw=.8)
                        ax.set_xticks([0,1],['OFF','ON'])
                    else:
                        for x,g in enumerate(['HC','OFF','ON']):
                            y=table.values[idx[g],k];y=y[np.isfinite(y)]
                            if len(y)>1 and np.ptp(y)>0:
                                parts=ax.violinplot(y,positions=[x],showmeans=True,showextrema=False)
                                for body in parts['bodies']: body.set_facecolor(colors[x])
                            ax.scatter(np.full(len(y),x),y,s=9,color=colors[x])
                        ax.set_xticks([0,1,2],['HC','OFF','ON'])
                    ax.set_title(table.names[k]);ax.set_ylabel(table.units[k])
                for ax in list(axes.flat)[len(subset):]: ax.set_visible(False)
                fig.tight_layout()
                # 使用数字编号作为文件名，避免直接将用户提供的族名用作路径。
                stem=kind if len(families)==1 and len(selected)<=16 else f'{kind}_family{families.index(family)+1}_{start//16+1}'
                for ext in ('png','pdf'):fig.savefig(output/(stem+'.'+ext),dpi=160)
                plt.close(fig)


def write_table(path,rows):
    with Path(path).open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def save_results(table, output, *, input_path=None, features_path=None):
    """保存统计表、受试者指标表、方法元数据和 PNG/PDF 图；拒绝覆盖已有目录。"""
    rows = compare_features(table)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    write_table(output / 'statistics.csv', rows)
    values = [dict(subject=scan['subject'], condition=scan['condition'],
                   **dict(zip(table.names, row))) for scan, row in zip(table.scans, table.values)]
    write_table(output / 'subject_features.csv', values)
    metadata = asdict(table)
    metadata.pop('values')
    metadata.update(input=str(Path(input_path).resolve()) if input_path is not None else None,
        features=str(Path(features_path).resolve()) if features_path is not None else None,
        correction='BH within each contrast and declared family; undefined tests retain slots',
        confidence_intervals='95%, unadjusted', missing='pairwise for ON/OFF, featurewise for HC comparisons')
    (output / 'metadata.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding='utf-8')
    plot_features(table, output)
    return output
