# https://www.mdpi.com/2673-6284/12/3/56
# https://www.nature.com/articles/s41598-018-24298-y

import itertools
import warnings
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact, f_oneway, tukey_hsd, kruskal, mannwhitneyu, ttest_ind, entropy, pearsonr#, false_discovery_control
from scipy.stats.contingency import odds_ratio
import scikit_posthocs as sp

def mannwhitney_test(df, only_pvalue:bool = True, pvalue_filter_limit_wilcox:float = 0.05, return_only_significant:bool = True):
    """
    Perform the Mann-Whitney U rank test on two independent samples between two conditions.

    Parameters
    ----------
    df: df
        pandas DataFrame of editing frequency. Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions.

    Returns
    -------
    DataFrame
        returns p-values for Mann-Whitney U rank test.
    """
    if df.iloc[:,:-2].empty:
        print('The input is an empty DataFrame. In this case the function returns empty list []')
        return False

    res = []
    pos = []
    index_comb = []
    cols_mw = list(itertools.combinations( df.iloc[:,-1].unique(), 2) )
    conditions = df.iloc[:,-1].unique()

    for c in df.columns[:-2]:
        data_by_condition = []
        for cond in conditions:
            data_by_condition.append(df[df.iloc[:,-1]==cond][c].values)
        aux_mw = mannwhitneyu(*data_by_condition)
        if (return_only_significant) and (aux_mw.pvalue <= pvalue_filter_limit_wilcox):
            if only_pvalue:
                res.append(aux_mw.pvalue)
                pos.append(c)
            else:
                res.append(aux_mw)
                pos.append(c)
        elif not return_only_significant:
            if only_pvalue:
                res.append(aux_mw.pvalue)
                pos.append(c)
            else:
                res.append(aux_mw)
                pos.append(c)

    df_pv = pd.DataFrame([pos,res], index=['coord','pvalue']).T.set_index('coord').T#.pivot(index='tests',columns='coord', values='p_value')
    df_pv[df.iloc[:,-2].name] = df.iloc[:,-2].unique()
    df_pv[df.iloc[:,-1].name] = '_vs_'.join(conditions)

    return df_pv



def t_student_test(df, only_pvalue:bool = True, pvalue_filter_limit_ttest:float = 0.05, return_only_significant:bool = True):
    """
    Perform the t Student test on two independent samples between two conditions.

    Parameters
    ----------
    df: df
        pandas DataFrame of editing frequency. Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions.

    Returns
    -------
    DataFrame
        returns p-values for t Student test.
    """
    if df.iloc[:,:-2].empty:
        print('The input is an empty DataFrame. In this case the function returns empty list []')
        return False

    res = []
    pos = []
    index_comb = []
    cols_mw = list(itertools.combinations( df.iloc[:,-1].unique(), 2) )
    conditions = df.iloc[:,-1].unique()

    for c in df.columns[:-2]:
        data_by_condition = []
        for cond in conditions:
            data_by_condition.append(df[df.iloc[:,-1]==cond][c].values)
        aux_ts = ttest_ind(*data_by_condition)
        if (return_only_significant) and (aux_ts.pvalue <= pvalue_filter_limit_ttest):
            if only_pvalue:
                res.append(aux_ts.pvalue)
                pos.append(c)
            else:
                res.append(aux_ts)
                pos.append(c)
        elif not return_only_significant:
            if only_pvalue:
                res.append(aux_ts.pvalue)
                pos.append(c)
            else:
                res.append(aux_ts)
                pos.append(c)

    df_pv = pd.DataFrame([pos,res], index=['coord','pvalue']).T.set_index('coord').T#.pivot(index='tests',columns='coord', values='p_value')
    df_pv[df.iloc[:,-2].name] = df.iloc[:,-2].unique()
    df_pv[df.iloc[:,-1].name] = '_vs_'.join(conditions)

    return df_pv

def anova_tukey_test(df, only_pvalue:bool = True, pvalue_filter_limit_anova:float = 0.05, pvalue_filter_limit_tukey:float = 0.05, return_only_significant:bool = True):
    """
    # Need to test in more than two conditions and to return a dataframe
    Anova with post-hoc test for more than two conditions.

    Parameters
    ----------
    df: df
        pandas DataFrame of editing frequency. Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions.

    Returns
    -------
    DataFrame
        returns p-values for Anova post-hoc test.
    """
    if df.iloc[:,:-2].empty:
        print('The input is an empty DataFrame. In this case the function returns empty list []')
        return False
    res = []
    pos = []
    index_comb = []
    cols_tukey = list(itertools.combinations( df.iloc[:,-1].unique(), 2) )
    conditions = df.iloc[:,-1].unique()
    for c in df.columns[:-2]:
        data_by_condition = []
        for cond in conditions:
            data_by_condition.append(df[df.iloc[:,-1]==cond][c].values)
        aux_aov = f_oneway(*data_by_condition)
        if aux_aov.pvalue<=pvalue_filter_limit_anova:
            aux_tukey = tukey_hsd(*data_by_condition)
            if return_only_significant and ((aux_tukey.pvalue <= pvalue_filter_limit_tukey).astype(int).sum() > 0):
                for cond_comb in itertools.combinations(range(len(conditions)),2):
                    if only_pvalue:
                        res.append(aux_tukey.pvalue[cond_comb])
                        pos.append(c)
                        index_comb.append(cond_comb)
                    else:
                        res.append(aux_tukey)
                        pos.append(c)
                        index_comb.append(cond_comb)
            elif not return_only_significant:
                for cond_comb in itertools.combinations(range(len(conditions)),2):
                    if only_pvalue:
                        res.append(aux_tukey.pvalue[cond_comb])
                        pos.append(c)
                        index_comb.append(cond_comb)
                    else:
                        res.append(aux_tukey)
                        pos.append(c)
                        index_comb.append(cond_comb)

    comb = [(conditions[c[0]], conditions[c[1]]) for c in index_comb]
    df_pv = pd.DataFrame([comb,pos,res], index=['tests','coord','pvalue']).T.pivot(index='tests',columns='coord', values='pvalue')
    return df_pv

def kruskal_dunn_test(df, only_pvalue:bool = True, pvalue_filter_limit_kruskal:float = 0.05, pvalue_filter_limit_dunn:float = 0.05, return_only_significant:bool = True, p_adjust:str = None):
    """
    # Still NEED to test in more than two conditions and to return a dataframe
    Kruskal-Wallis with post-hoc test for more than two conditions.

    Parameters
    ----------
    df: df
        pandas DataFrame of editing frequency. Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions.

    Returns
    -------
    DataFrame
        returns p-values for Anova post-hoc test.
    """
    if df.iloc[:,:-2].empty:
        print('The input is an empty DataFrame. In this case the function returns empty list []')
        return False
    res = []
    pos = []
    index_comb = []
    cols_tukey = list(itertools.combinations( df.iloc[:,-1].unique(), 2) )
    conditions = df.iloc[:,-1].unique()
    group_col = df.columns[-1]
    for c in df.columns[:-2]:
        data_by_condition = []
        for cond in conditions:
            data_by_condition.append(df[df.iloc[:,-1]==cond][c].values)
        aux_aov = kruskal(*data_by_condition)
        if aux_aov.pvalue<=pvalue_filter_limit_kruskal:
            # aux_tukey = tukey_hsd(*data_by_condition)
            aux_dunn = sp.posthoc_dunn(a=df, val_col=c, group_col=group_col, p_adjust=p_adjust)
            if return_only_significant and ( (aux_dunn<=pvalue_filter_limit_dunn).astype(int).sum().sum()>0 ):
                for cond_comb in itertools.combinations(conditions, 2):
                    res.append(aux_dunn.loc[cond_comb])
                    pos.append(c)
                    index_comb.append(cond_comb)
            elif not return_only_significant:
                for cond_comb in itertools.combinations(conditions, 2):
                    res.append(aux_dunn.loc[cond_comb])
                    pos.append(c)
                    index_comb.append(cond_comb)

    # comb = [(conditions[c[0]], conditions[c[1]]) for c in index_comb]
    df_pv = pd.DataFrame([index_comb,pos,res], index=['tests','coord','pvalue']).T.pivot(index='tests',columns='coord', values='pvalue')
    return df_pv

def chi2_test(df_a, df_g, only_pvalue=True, return_only_significant=True, pvalue_filter_limit=0.05):
    aux_pv = []
    for c in df_a.columns[:-2]:
        if only_pvalue:
            aux_pv.append(chi2_contingency([df_a[c].values, df_g[c].values], lambda_="log-likelihood").pvalue)
        else:
            aux_pv.append(chi2_contingency([df_a[c].values, df_g[c].values], lambda_="log-likelihood"))

    if only_pvalue:
        chi2_res = pd.DataFrame(aux_pv, index=df_a.columns[:-2], columns=['pvalue'])
    else:
        chi2_res = pd.DataFrame(aux_pv, index=df_a.columns[:-2], columns=['statistic','pvalue','dof','expected_freq'])

    if return_only_significant:
        chi2_res = chi2_res[chi2_res.pvalue<=pvalue_filter_limit].T

    chi2_res[df_a.columns[-2]] = df_a.iloc[0,-2]
    chi2_res[df_a.columns[-1]] = '_vs_'.join(df_a.iloc[:,-1].unique())
    return chi2_res

def fisher_test(df_a, df_g, only_pvalue=True, return_only_significant=True, pvalue_filter_limit=0.05):
    aux_pv = []
    for c in df_a.columns[:-2]:
        if only_pvalue:
            aux_pv.append(fisher_exact([df_a[c].values, df_g[c].values], alternative='two-sided').pvalue)
        else:
            aux_pv.append(fisher_exact([df_a[c].values, df_g[c].values], alternative='two-sided'))

    if only_pvalue:
        fisher_res = pd.DataFrame(aux_pv, index=df_a.columns[:-2], columns=['pvalue'])
    else:
        fisher_res = pd.DataFrame(aux_pv, index=df_a.columns[:-2], columns=['statistic','pvalue'])

    if return_only_significant:
        fisher_res = fisher_res[fisher_res.pvalue<=pvalue_filter_limit].T
    else:
        fisher_res = fisher_res.T

    fisher_res[df_a.columns[-2]] = df_a.iloc[0,-2]
    fisher_res[df_a.columns[-1]] = '_vs_'.join(df_a.iloc[:,-1].unique())

    return fisher_res

def odds_r(df_a, df_g):
    aux_or = []
    for c in df_a.columns[:-2]:
        aux_or.append(odds_ratio([df_a[c].astype(int).values, df_g[c].astype(int).values], kind='conditional').statistic)
    odds_res = pd.DataFrame(aux_or, index=df_a.columns[:-2], columns=['odds_ratio']).T

    odds_res[df_a.columns[-2]] = df_a.iloc[0,-2]
    odds_res[df_a.columns[-1]] = '_vs_'.join(df_a.iloc[:,-1].unique())

    return odds_res

def entropy_calculation(df_a, df_g):
    a = df_a.copy()
    g = df_g.copy()
    condition_name = a.columns[-1]
    conditions = a[condition_name].unique()
    cols = a.columns[:-2]
    a['base_aux'] = 'a'
    g['base_aux'] = 'g'

    aux = pd.concat([a, g])
    del(a)
    del(g)
    aux_dict = {cond:[] for cond in conditions}

    for cond in conditions:
        for col in aux.columns[:-3]:
            aux_dict[cond].append(entropy(aux[(aux[condition_name]==cond)][col].values/aux[(aux[condition_name]==cond)][col].values.sum()))
    aux = pd.DataFrame(aux_dict).T
    aux.columns = cols
    return aux

def conditon_pearson_corr(df, return_only_significant:bool = True, pvalue_filter_limit:float = 0.05):
    aux = df.copy()

    condition_name = aux.columns[-1]
    conditions = aux[condition_name].unique()
    P_list = []
    pv_list = []

    for col in aux.columns[:-2]:
        P_list.append(pearsonr(aux[aux[condition_name]==conditions[0]][col], aux[aux[condition_name]==conditions[1]][col]).statistic)
        pv_list.append(pearsonr(aux[aux[condition_name]==conditions[0]][col], aux[aux[condition_name]==conditions[1]][col]).pvalue)
    aux = pd.DataFrame([P_list,pv_list], columns=aux.columns[:-2], index=['R','pvalue']).T.reset_index()

    if return_only_significant:
        aux = aux[aux.pvalue <= pvalue_filter_limit]

    aux.columns = ['Positions', 'R', 'pvalue']

    return aux

def gene_pearson_corr(df, exprs:list, return_only_significant:bool = True, pvalue_filter_limit:float = 0.05):
    aux = df.copy()

    P_list = []
    pv_list = []

    for col in aux.columns[:-2]:
        P_list.append(pearsonr(aux[col].values, exprs).statistic)
        pv_list.append(pearsonr(aux[col].values, exprs).pvalue)
    aux = pd.DataFrame([P_list,pv_list], columns=aux.columns[:-2], index=['R','pvalue']).T.reset_index()

    if return_only_significant:
        aux = aux[aux.pvalue <= pvalue_filter_limit]
    aux.columns = ['Positions', 'R', 'pvalue']
    return aux
