# https://www.mdpi.com/2673-6284/12/3/56
# https://www.nature.com/articles/s41598-018-24298-y

import itertools
import warnings
# import pandas as pd
# import numpy as np
# from scipy.stats import chi2_contingency, fisher_exact, f_oneway, tukey_hsd, kruskal
# from scipy.stats.contingency import odds_ratio
# import scikit_posthocs as sp

def merge_files_one_region(meta):
    """
    Merge all RES files for the same and UNIQUE region (output of REDItools2) in three pandas DataFrames of frequency or count per position.
    The first DataFrame is the frequency of editing (A-G or T-C). The Second DataFrame is the count of A (or T) per position. And the last one is the count of G (or C) per postion.

    Parameters
    ----------
    meta: df
        A pandas DataFrame with metadata information for UNIQUE region. The first four columns are mandatory
            First: Full path file names of REDItools2 results tables
            Second: Samples names
            Third: region (gene symbol)
            Fourth: conditions

    Returns
    -------
    tuple
        a tuple of three pandas DataFrames (df, df_a, df_g).
        df: frequency of editing
        df_a: counts of A or T
        df_g: counts of G or C
    """
    FILES = meta.iloc[:,0].unique()
    # region = meta.iloc[:,2].values[0]
    df_list_a = []
    df_list_g = []
    samples = []
    for file in FILES:
        df = pd.read_csv( file, sep = '\t' )
        samples.append(meta[meta.iloc[:,0] == file].iloc[:,1].values[0])
        if df.empty:
            df_list_a.append( pd.DataFrame({'Position':[-1],
                                            samples[-1]:[np.nan]}).set_index('Position') )
            df_list_g.append( pd.DataFrame({'Position':[-1],
                                            samples[-1]:[np.nan]}).set_index('Position') )
        else:
            df = df[df['Coverage-q30']>=10]
            df = df[df.Reference.isin(['A','T'])]
            df = df[(df.AllSubs.str.contains('AG')) | (df.AllSubs.str.contains('TC'))]
            if df.empty:
                df_list_a.append( pd.DataFrame({'Position':[-1], samples[-1]:[np.nan]}).set_index('Position') )
                df_list_g.append( pd.DataFrame({'Position':[-1], samples[-1]:[np.nan]}).set_index('Position') )
            else:
                df[['count_A','count_C','count_G','count_T']] = df['BaseCount[A,C,G,T]'].str.replace('[','',regex=False).str.replace(']','', regex=False).str.split(',',expand=True).apply(pd.to_numeric)

                df[samples[-1]] = np.where((df['Reference']=='A'), df.count_A, df.count_T)
                df_list_a.append(df[['Position',samples[-1]]].set_index('Position'))

                df[samples[-1]] = np.where((df['Reference']=='A'), df.count_G, df.count_C)
                df_list_g.append(df[['Position',samples[-1]]].set_index('Position'))
    # Analysis by gene
    df_a = pd.concat(df_list_a, axis=1)
    df_g = pd.concat(df_list_g, axis=1)

    if (-1 in df_a.index):
        df_a.drop(-1, axis=0, inplace=True)
        df_g.drop(-1, axis=0, inplace=True)

    if df_a.empty:
        # print('Region', region , 'has not filterd edited site.')
        # return df_a.T
        warnings.warn("Sorry, all dataset are empty or didn't has A-to-G editing.")
        return df.T, df_a.T, df_g.T
    else:
        df = pd.DataFrame(100*df_g.values / (df_a.values+df_g.values), columns=samples, index=df_a.index).T
        df = df.merge(meta.iloc[:,[1,2,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)

        df_a = df_a.T.merge(meta.iloc[:,[1,2,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)
        df_g = df_g.T.merge(meta.iloc[:,[1,2,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)
        return df, df_a, df_g

def merge_files_all_regions(meta):
    """
    Merge all RES files for ALL regions (output of REDItools2) in three pandas DataFrames of frequency or count per position.
    The first DataFrame is the frequency of editing (A-G or T-C). The Second DataFrame is the count of A (or T) per position. And the last one is the count of G (or C) per postion.

    Parameters
    ----------
    meta: df
        A pandas DataFrame with metadata information for ALL regions of interest. The first four columns are mandatory
            First: Full path file names of REDItools2 results tables
            Second: Samples names
            Third: regions (genes symbol)
            Fourth: conditions

    Returns
    -------
    tuple
        a tuple of three pandas DataFrames (df, df_a, df_g).
        df: frequency of editing
        df_a: counts of A or T
        df_g: counts of G or C
        region_list: list of non-empty regions
    """
    df_list = []
    df_a_list = []
    df_g_list = []
    region_list = dict()
    df, df_a, df_g = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    for r in meta.iloc[:,2].unique():
        # print(r+'+ ', end=' ')
        try:
            df, df_a, df_g = merge_files_one_region(meta[meta.iloc[:,2]==r])
            # print(df+'- ', end=' ')
            df_list.append(df.iloc[:,:-2])
            df_a_list.append(df_a.iloc[:,:-2])
            df_g_list.append(df_g.iloc[:,:-2])
            region_list[r] = df.columns[:-2]

        except:
            pass
    # return df_list, df_a_list, df_g_list, region_list

    if not df_list:
        warnings.warn("All dataset for region "+r+"are empty or didn't has A-to-G aditing.")
    elif not df_a_list:
        warnings.warn("All dataset for region "+r+"are empty or didn't has A-to-G aditing.")
    elif not df_g_list:
        warnings.warn("All dataset for region "+r+"are empty or didn't has A-to-G aditing.")

    try:
        df = pd.concat(df_list, axis=1).merge(meta.iloc[:,[1,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)
        df_a = pd.concat(df_a_list, axis=1).merge(meta.iloc[:,[1,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)
        df_g = pd.concat(df_g_list, axis=1).merge(meta.iloc[:,[1,3]].drop_duplicates().set_index(meta.columns[1]), left_index=True, right_index=True)

        v = ['several' for r in range(df.shape[0])]
        i = len(df.columns) - 1
        df.insert( i, 'region', v)
        df_a.insert( i, 'region', v)
        df_g.insert( i, 'region', v)

        return df, df_a, df_g, region_list
    except:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), region_list

def filter_snps(df, vcf_file=''):
    pass

def filter_positions(df, nan_filter=True, nan_filter_limit=0,
                     zero_filter=True, zero_filter_limit=0,
                     hundred_filter=True, hundred_filter_limit=0):
    """
    Filter positions limiting the quantity of samples with nan values, zero editing and 100% of editing across samples.

    Parameters
    ----------
    meta: df
        Frequency editing DataFrame of all samples analyzed. The DataFrame must be like merge_files output. The last two columns must be region and conditions.
    nan_filter: bool
        Must be True to filter quantity of nan values across samples in each contition
    nan_filter_limit: int
        Number max of samples with nan in one position in each contition
    zero_filter: bool
        Must be True to filter quantity of zero editing values across samples in each contition
    zero_filter_limit: int
        Number max of samples with zero frequency in one position in each contition
    hundred_filter: bool
        Must be True to filter quantity of 100% of editing values across samples in each contition
    hundred_filter_limit: int
        Number max of samples with 100% frequency in one position in each contition

    Returns
    -------
    df
        DataFrame without filtered positions
    """
    if df.iloc[:,:-2].empty:
        # print('The input is an empty DataFrame. In this case the function returns the input')
        return df
    array_filter = df.iloc[:,range(df.shape[1]-1)].columns # minus 1 to exclude type column
    for cond in df.iloc[:,-1].unique():
        # nan filter
        if nan_filter:
            array_filter = array_filter[ np.where( df.loc[df.iloc[:,-1]==cond, array_filter ].isna().sum(axis=0)<=nan_filter_limit )[0]]
        # zero filter
        if zero_filter:
            array_filter = array_filter[ np.where( (df.loc[df.iloc[:,-1]==cond, array_filter ]==0).sum(axis=0)<=zero_filter_limit )[0]]
        # hundred filter
        if hundred_filter:
            array_filter = array_filter[ np.where( (df.loc[df.iloc[:,-1]==cond, array_filter ]==100).sum(axis=0)<=hundred_filter_limit )[0]]
        # G-test filter power_divergence http://www.biostathandbook.com/repgtestgof.html
        # if g_test_filter:
        #     array_filter = array_filter[ np.where( (df.loc[df.iloc[:,-1]==cond, array_filter ]==100).sum(axis=0)<=hundred_filter_limit )[0]]

    # var filter
    # array_filter = array_filter[ np.where((df.loc[df.iloc[:,-1]==cond, array_filter ].var()>=.01))[0]]

    array_filter = array_filter.tolist()
    array_filter.append( df.columns[-1] )
    return df[array_filter]

def independency_gtest(df_a, df_g, only_pvalue=True):
    """
    Perform independecy G-test to verify if there's deviation from the expected proportions and significant variation among the replicates.

    Parameters
    ----------
    df_a: df
        pandas DataFrame of Adenine counts (A). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    df_g: df
        pandas DataFrame of Guanine (Inosine) counts (G). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    only_pvalue: bool
        If True, the function will return two DataFrames only with p-values. Otherwise will return the attributes of Chi2ContingencyResult object (https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.chi2_contingency.html).

    Returns
    -------
    df
        DataFrame with scipy.stats.chi2_contingency results.
    """
    p_values = {c:[] for c in df_a.iloc[:,-1].unique()}
    if (df_a.iloc[:,:-2].empty) or (df_g.iloc[:,:-2].empty):
        # print('The input is an empty DataFrame. In this case the function returns the input')
        return df_a
    array_filter = df_a.iloc[:,range(df_a.shape[1]-2)].columns
    for cond in p_values.keys():
        for arr in array_filter:
            a = df_a.loc[df_a.iloc[:,-1]==cond, arr ].values
            g = df_g.loc[df_g.iloc[:,-1]==cond, arr ].values
            if only_pvalue:
                if (a.sum() <= 0) or (g.sum() <= 0): # only zeros raise a Error in chi2_contingency
                    p_values[cond].append(np.nan)
                else:
                    p_values[cond].append(chi2_contingency(np.array([a,g]), lambda_="log-likelihood").pvalue)
            else:
                if (a.sum() <= 0) or (g.sum() <= 0): # only zeros raise a Error in chi2_contingency
                    p_values[cond].append(np.nan)
                else:
                    p_values[cond].append(chi2_contingency(np.array([a,g]), lambda_="log-likelihood"))
    aux = pd.DataFrame(p_values, index=array_filter).T.dropna(axis=1)
    return aux

def filter_gtest(df_a, df_g, pvalue_filter_limit=0.05, gtest_filter_limit=0, bh_correction=False):
    """
    Filter positions limiting the quantity of samples with significant p-values for independecy G-test.

    Parameters
    ----------
    df_a: df
        pandas DataFrame of Adenine counts (A). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    df_g: df
        pandas DataFrame of Guanine (Inosine) counts (G). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    pvalue_filter_limit: float
        Minimum p-value threshold to consider the chi2_test significant.
    gtest_filter_limit: int
        Number max of conditions with significant p-value frequency in one position in each contition
    bh_correction: bool
        If True, all p-values will be corrected by false_discovery_control

    Returns
    -------
    tuple
        returns inputs df_a, df_g without filtered coordinates.
    """
    aux = independency_gtest(df_a, df_g, only_pvalue=True)
    if bh_correction:
        aux.iloc[:,:] = false_discovery_control(aux, axis=1)

    array_filter = aux.columns
    array_filter = array_filter[ np.where((aux <= pvalue_filter_limit).astype(int).sum(axis=0)<=gtest_filter_limit)[0] ]

    array_filter = array_filter.tolist()
    array_filter.append( df_a.columns[-2] )
    array_filter.append( df_a.columns[-1] )
    return df_a[array_filter], df_g[array_filter]

def pool_positions(df_a, df_g, pvalue_filter_limit=0.05, gtest_filter_limit=0, bh_correction=False):
    """
    Pool the counts per positions and per consitions of Adenine (A) and Guanine (G) Dataframes. To improve the Fisher exact test or Chi2 test between conditions.

    Parameters
    ----------
    df_a: df
        pandas DataFrame of Adenine counts (A). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    df_g: df
        pandas DataFrame of Guanine (Inosine) counts (G). Rows are samples and columns are coordinates. The DataFrame must be like merge_files output. The last two columns must be region and conditions
    pvalue_filter_limit: float
        Minimum p-value threshold to consider the chi2_test significant.
    gtest_filter_limit: int
        Number max of conditions with significant p-value frequency in one position in each contition
    bh_correction: bool
        If True, all p-values will be corrected by false_discovery_control

    Returns
    -------
    tuple
        returns inputs df_a, df_g without filtered coordinates.
    """
    aux_a, aux_g = filter_gtest(df_a, df_g, pvalue_filter_limit, gtest_filter_limit, bh_correction)

    aux_a = aux_a.groupby(df_a.columns[-1]).sum()
    aux_a[df_a.columns[-2]] = df_a.iloc[0,-2]
    aux_a[df_a.columns[-1]] = aux_a.index

    aux_g = aux_g.groupby(df_g.columns[-1]).sum()
    aux_g[df_g.columns[-2]] = df_g.iloc[0,-2]
    aux_g[df_g.columns[-1]] = aux_g.index

    return aux_a, aux_g

def wilcox_test(df, only_pvalue:bool = True, pvalue_filter_limit_wilcox:float = 0.05, return_only_significant:bool = True):
    pass

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
    df_pv = pd.DataFrame([comb,pos,res], index=['tests','coord','p_value']).T.pivot(index='tests',columns='coord', values='p_value')
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
    df_pv = pd.DataFrame([index_comb,pos,res], index=['tests','coord','p_value']).T.pivot(index='tests',columns='coord', values='p_value')
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
# def chi2_test(df_a, df_b):
#     power_divergence(table) para todas posiçoes e tipos para avaliar a divergencia
    # return []

#                 res.append(pg.anova(dv=c, between=df.columns[-1], data=df, detailed=False))#.to_csv('A2I_anova.csv')
#                 res[-1]['Position'] = c

#                 aux_tukey = aux_tukey.loc[df.iloc[:,-1].unique(), df.iloc[:,-1].unique()] # just to confirm that the order will be always the same
#                 np.fill_diagonal(aux_tukey.values, 0) # to transform to squareform, the diagonal must be zero
#                 res[-1].loc[:,cols_tukey] = squareform(aux_tukey)
    #
#     if res:
#         res = pd.concat(res).set_index('Position')
#         # stat_res.append(res.merge(df[res.index].T, left_index=True, right_index=True))
#         # stat_res[-1]['gene'] = g

#         return aux_tukey
#     else:
#         return False
# sep = '\t'
# def merge_files(meta, sep='\t'):
#     REGIONS = meta.iloc[:,2].unique()[:1]


# def merge_utr(res, utr):
#     utr_list = []
#     for i in range(res.shape[0]):
#         p = res.iloc[i, 7]
#         r = res.iloc[i, 8]
#         if not utr[(utr.gene_name==r) & (utr.utr_max >= p) & (utr.utr_min <= p)].empty:
#             utr_list.append(utr[(utr.gene_name==r) & (utr.utr_max>=p) & (utr.utr_min<=p)].utr_type.values[0])
#         else:
#             utr_list.append('non_utr')
#     res['utr_type'] = utr_list
#     return res
