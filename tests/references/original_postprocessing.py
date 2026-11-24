# ruff: noqa: RUF059
# 保留上游接收但未使用的 p_value 变量，不改参考函数计算语句。
"""原作者 Python 后处理参考；仅封装输入输出，不改变统计公式。

出处： https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/group_comparison.py
功能分离度： https://github.com/chenfei-ye/WEIDA-DBS/blob/2e20b0bc956e6eb0eb8e06e2aaf00a93b014e72c/calcu_FS.py
"""
import math

import numpy as np
import scipy.stats
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests


def get_fdr(p_list):
    fdr = multipletests(p_list, method="fdr_bh", alpha=0.05)[1]

    return fdr

def calculate_functional_segregation_strength(ts, modu1, modu2):
    fc1=[]
    fc2=[]
    fc12=[]
    for i in range(len(modu1)):
        for j in range(i+1,len(modu1)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu1[i]],ts.iloc[:,modu1[j]])
            fc1.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
    for i in range(len(modu2)):
        for j in range(i+1,len(modu2)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu2[i]],ts.iloc[:,modu2[j]])
            fc2.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
    for i in range(len(modu1)):
        for j in range(len(modu2)):
            corr_coeff, p_value = pearsonr(ts.iloc[:,modu1[i]],ts.iloc[:,modu2[j]])
            fc12.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))

    within_module_fc = fc1 + fc2
    across_module_fc = fc12

    average_within_module_fc = np.mean(within_module_fc)

    average_across_module_fc = np.mean(across_module_fc)

    functional_segregation_strength = average_within_module_fc - average_across_module_fc

    return average_within_module_fc, average_across_module_fc, functional_segregation_strength

def calculate_nodal_functional_segregation_strength(ts, modu1, modu2, num_parcel):

    node_fs = [None] * num_parcel

    for i in range(num_parcel):
        fc_in = []
        fc_ac = []
        if i in modu1:
            for j in modu1:
                if j != i:
                    corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                    fc_in.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
            for j in modu2:
                corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                fc_ac.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
        if i in modu2:
            for j in modu1:
                corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                fc_ac.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))
            for j in modu2:
                if j != i:
                    corr_coeff, p_value = pearsonr(ts.iloc[:,i],ts.iloc[:,j])
                    fc_in.append(0.5*math.log((1+corr_coeff)/(1-corr_coeff)))

        node_fs[i] = np.mean(fc_in) - np.mean(fc_ac)

    return node_fs

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
