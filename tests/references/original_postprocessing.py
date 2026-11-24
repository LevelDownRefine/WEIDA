"""原作者 Python 后处理参考；仅封装输入输出，不改变统计公式。

出处： https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py
"""

import numpy as np
import scipy.stats
from statsmodels.stats.multitest import multipletests


def get_fdr(p_list):
    fdr = multipletests(p_list, method="fdr_bh", alpha=0.05)[1]

    return fdr

def original_group_comparison(hc, off, on):
    """原脚本 DBS 三组检验段；仅把 CSV 切片与固定 4 列改为函数输入。"""
    t_values, p_values = [], []
    for i in range(hc.shape[1]):
        c0 = list(hc[:, i])
        c1 = list(off[:, i])
        c2 = list(on[:, i])
        t0, pval0 = scipy.stats.ttest_ind(c1, c0)
        t1, pval1 = scipy.stats.ttest_rel(c2, c1)
        t2, pval2 = scipy.stats.ttest_ind(c2, c0)
        t_values.append([t0, t1, t2])
        p_values.append([pval0, pval1, pval2])
    t_values, p_values = np.array(t_values).T, np.array(p_values).T
    return t_values, p_values, np.array([get_fdr(p) for p in p_values])
