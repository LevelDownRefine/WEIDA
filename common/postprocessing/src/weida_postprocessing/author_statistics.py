"""作者 DBS 组间统计的正式实现，仅拆分步骤与参数。

出处：https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py
"""
import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests

GROUPS = ('HC', 'OFF', 'ON')
CONTRASTS = ('OFF-HC', 'ON-OFF', 'ON-HC')


def describe_groups(hc, off, on):
    """原第 82–87 行：均值、总体标准差；返回组 × 指标。"""
    means, stds = [], []
    for group in (hc, off, on):
        means.append([np.mean(list(group[:, i])) for i in range(group.shape[1])])
        stds.append([np.std(list(group[:, i])) for i in range(group.shape[1])])
    return np.array(means), np.array(stds)


def compare_groups(hc, off, on):
    """原第 89–91 行：三组 t 检验；返回对比 × 指标的 t、p。"""
    t_values, p_values = [], []
    for i in range(hc.shape[1]):
        c0, c1, c2 = list(hc[:, i]), list(off[:, i]), list(on[:, i])
        t0, p0 = scipy.stats.ttest_ind(c1, c0)
        t1, p1 = scipy.stats.ttest_rel(c2, c1)
        t2, p2 = scipy.stats.ttest_ind(c2, c0)
        t_values.append([t0, t1, t2])
        p_values.append([p0, p1, p2])
    return np.array(t_values).T, np.array(p_values).T


def get_fdr(p_list):
    """原第 6–9 行：对一个对比的指标做 BH 校正。"""
    return multipletests(p_list, method="fdr_bh", alpha=0.05)[1]


def group_comparison(hc, off, on):
    """串联检验和原第 115–117 行的逐对比校正，返回 t、p、q。"""
    t_values, p_values = compare_groups(hc, off, on)
    return t_values, p_values, np.array([get_fdr(p) for p in p_values])
