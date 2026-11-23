"""为 WEiDA 输出生成描述性图表，不进行显著性推断。"""
import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import ListedColormap


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    target=args.output or args.input/'figures'; target.mkdir(parents=True,exist_ok=True)
    available={f.name for f in font_manager.fontManager.ttflist}
    font=next((f for f in ['Microsoft YaHei','SimHei','Noto Sans CJK SC'] if f in available),'DejaVu Sans')
    plt.rcParams.update({'font.family':font,'axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    r=np.load(args.input/'results.npz'); meta=json.loads((args.input/'metadata.json').read_text(encoding='utf-8'))
    groups=['HC','OFF','ON']; names=['健康对照 HC','STN 刺激关闭 OFF','STN 刺激开启 ON']; colors=['#547aa5','#d98743','#479c8e']
    state_colors=['#5569ac','#d59a39','#469d8f','#b45e88']; conditions=np.array([x['condition'] for x in meta['scans']])
    counts=[int(np.sum(conditions==g)) for g in groups]
    note='探索性结果：统一特征向量符号约定；尚未验证 MATLAB 等价性。图中不表示组间差异具有统计显著性。'
    def save(fig,name):
        fig.savefig(target/(name+'.png'),dpi=180,bbox_inches='tight')
        fig.savefig(target/(name+'.pdf'),bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    fig.suptitle('WEiDA 四状态分析｜65 次扫描 · 78 个皮层脑区',fontsize=19)
    for ax,key,title,ylabel in [(axes[0,0],'occupancy','各组状态占比','占比（%）'),(axes[0,1],'lifetime','各组平均持续时间','秒')]:
        for gi,g in enumerate(groups):
            data=r[key][conditions==g]*(100 if key=='occupancy' else 1)
            avg=np.nanmean(data,axis=0); n=np.sum(np.isfinite(data),axis=0)
            sem=np.nanstd(data,axis=0,ddof=1)/np.sqrt(n)
            ax.bar(np.arange(4)+(gi-1)*.24,avg,width=.23,color=colors[gi],label=f'{names[gi]}（n={counts[gi]}）',yerr=sem,capsize=3)
        ax.set(xticks=range(4),xticklabels=['状态 1','状态 2','状态 3','状态 4'],ylabel=ylabel,title=title+'（扫描均值 ± 标准误）')
        ax.grid(axis='y',alpha=.15); ax.set_axisbelow(True)
    axes[0,0].legend(fontsize=8)
    ax=axes[1,0]
    for gi,g in enumerate(groups):
        data=r['occupancy'][conditions==g]*100
        for j in range(4):
            jitter=np.linspace(-.075,.075,len(data))
            ax.scatter(j+(gi-1)*.24+jitter,data[:,j],s=14,color=colors[gi],alpha=.65)
    ax.set(xticks=range(4),xticklabels=['状态 1','状态 2','状态 3','状态 4'],ylabel='占比（%）',title='每次扫描的状态占比（每点代表一次扫描）')
    ax.grid(axis='y',alpha=.15)
    ax=axes[1,1]; total=np.bincount(r['states'],minlength=4); pct=total/total.sum()*100
    bars=ax.bar(range(4),pct,color=state_colors)
    for bar,n,p in zip(bars,total,pct): ax.text(bar.get_x()+bar.get_width()/2,p+.5,f'{p:.1f}%\n{n:,} 帧',ha='center',fontsize=10)
    ax.set(xticks=range(4),xticklabels=['状态 1','状态 2','状态 3','状态 4'],ylabel='全部时间点占比（%）',title='全体状态分布（按频率排序）',ylim=(0,max(pct)*1.25))
    fig.supxlabel(note,fontsize=9)
    save(fig,'01_summary')
    fig,axes=plt.subplots(1,3,figsize=(13,4.7),layout='constrained')
    for ax,g,name in zip(axes,groups,names):
        m=r['group_'+g]; im=ax.imshow(m,vmin=0,vmax=1,cmap='Blues')
        for i in range(4):
            for j in range(4): ax.text(j,i,f'{m[i,j]:.3f}',ha='center',va='center',color='white' if m[i,j]>.55 else '#202832')
        ax.set(xticks=range(4),yticks=range(4),xticklabels=range(1,5),yticklabels=range(1,5),xlabel='下一状态',ylabel='当前状态',title=name)
    fig.colorbar(im,ax=axes,shrink=.7,label='转移概率')
    fig.suptitle('组平均转移概率｜包含状态保持（对角线）',fontsize=17)
    fig.supxlabel('先计算每次扫描的转移概率，再按有效扫描平均；状态编号在各组间共享。\n'+note,fontsize=9)
    save(fig,'02_transitions')
    fig,ax=plt.subplots(figsize=(18,5),layout='constrained'); limit=np.max(np.abs(r['centers']))
    im=ax.imshow(r['centers'],aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
    ax.set(yticks=range(4),yticklabels=['状态 1','状态 2','状态 3','状态 4'],xticks=range(78),xticklabels=meta['regions'],title='四个状态的加权特征向量中心｜AAL1ctx 78 脑区')
    ax.tick_params(axis='x',labelrotation=90,labelsize=6)
    fig.colorbar(im,ax=ax,shrink=.65,label='特征向量权重')
    fig.supxlabel('颜色是模型特征权重，不是 BOLD 激活强度；保留输入文件的脑区顺序。\n'+note,fontsize=9)
    save(fig,'03_state_patterns')
    lengths=r['lengths']; canvas=np.full((len(lengths),max(lengths)),np.nan); start=0
    for i,n in enumerate(lengths): canvas[i,:n]=r['states'][start:start+n]+1;start+=n
    fig,ax=plt.subplots(figsize=(13,9),layout='constrained'); cmap=ListedColormap(state_colors);cmap.set_bad('#ededed')
    im=ax.imshow(canvas,aspect='auto',interpolation='nearest',cmap=cmap,vmin=.5,vmax=4.5)
    for boundary in np.flatnonzero(conditions[1:]!=conditions[:-1])+.5: ax.axhline(boundary,color='white',linewidth=2)
    ax.set(yticks=range(len(lengths)),yticklabels=[x['subject'] for x in meta['scans']],xlabel=f'时间点索引（从 0 开始；TR={meta["tr"]} 秒）',title='全部扫描的状态序列｜每一行是一份独立扫描')
    ax.tick_params(axis='y',labelsize=6)
    fig.colorbar(im,ax=ax,ticks=[1,2,3,4],shrink=.4,label='状态')
    fig.supxlabel('灰色为该扫描不存在的时间点，不做填补。\n'+note,fontsize=9)
    save(fig,'04_state_sequences')
    print(target.resolve())

if __name__=='__main__': main()
