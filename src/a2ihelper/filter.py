# https://www.mdpi.com/2673-6284/12/3/56
# https://www.nature.com/articles/s41598-018-24298-y

import warnings
import requests
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency#, false_discovery_control
from scipy.stats.contingency import odds_ratio
import scikit_posthocs as sp

def merge_files_one_region(meta, coverage_q30:int = 10):
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

    coverage_q30: int
        An integer to define the minimum of reads with quality higher than in each position

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
        df['Position'] = df['Region'].str.replace('chr','') + '_' + df['Position'].astype(str)
        samples.append(meta[meta.iloc[:,0] == file].iloc[:,1].values[0])
        if df.empty:
            df_list_a.append( pd.DataFrame({'Position':[-1],
                                            samples[-1]:[np.nan]}).set_index('Position') )
            df_list_g.append( pd.DataFrame({'Position':[-1],
                                            samples[-1]:[np.nan]}).set_index('Position') )
        else:
            df = df[df['Coverage-q30']>=coverage_q30]
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

def merge_files_all_regions(meta, coverage_q30:int = 10):
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

    coverage_q30: int
        An integer to define the minimum of reads with quality higher than in each position

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
            df, df_a, df_g = merge_files_one_region(meta[meta.iloc[:,2]==r], coverage_q30=coverage_q30)
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

def call_snp_vep(coordinates:list = [], species:str = 'homo_sapiens', sub:str = 'G'):
    """
    Consult SNPs in VEP (ensembl tool).

    Parameters
    ----------
    coordinates: list
        List of strings to call VEP. The string must be a list of chr_coordinate (e.g: ['9_129401662', '1_6524705']).
    species: str
        String with the species that you are requesting. It is based on Ensembl names like homo_sapiens, mus_musculus, rattus_norvegicus, zebrafish. Other genomes can be searched here: https://rest.ensembl.org/documentation/info/species
    sub: str
        Substitution base, can be G (A>G) or C (T>C).

    Returns
    -------
    list
        List of chr_coordinate and respective rsID
    """

    if sub == 'G':
        ref = 'A'
    elif sub == 'C':
        ref = 'T'
    else:
        raise Exception("*Base not allowed* -> Sorry, only G or C are allowed as value for 'sub' parameter.")
    server = 'https://rest.ensembl.org'
    aux = []
    for coord in coordinates:
        chr, crd = coord.split('_')
        ext = '/vep/'+species+'/region/'+chr+':'+crd+"-"+crd+'/'+sub+'?'
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if r.ok:
            decoded = r.json()
            for d in decoded:
                if 'colocated_variants' in d.keys():
                    for d_ in d['colocated_variants']:
                        if d_['allele_string'].startswith(ref):
                            aux.append([chr+'_'+crd, d_['id']])
    if not aux:
        warnings.warn('*Returning empty list* -> Positions of genes were not found in the VEP annotation as SNP.')
    return aux

def filter_snps_vep(df, species:str = 'homo_sapiens', sub:str = 'G'):
    """
    Filter positions by consulting SNPs in VEP (ensembl tool).

    Parameters
    ----------
    df: DataFrame
        Frequency editing DataFrame of all samples analyzed. The DataFrame must be like merge_files output. The columns must be chr_coordinate (e.g: 9_129401662). The last two columns must be region and conditions.
    species: str
        String with the species that you are requesting. It is based on Ensembl names like homo_sapiens, mus_musculus, rattus_norvegicus, zebrafish. Other genomes can be searched here: https://rest.ensembl.org/documentation/info/species
    sub: str
        Substitution base, can be G (A>G) or C (T>C).

    Returns
    -------
    df_1
        DataFrame df without filtered positions

    df_2
        DataFrame with positions and rsID
    """

    aux = call_snp_vep(df.iloc[:,:-2].columns, species, sub)
    aux = pd.DataFrame(aux, columns=['Position','rsID'])
    return df.loc[:,~df.columns.isin(aux['Position'].values)], aux

def filter_positions(df, nan_filter=True, nan_filter_limit=0,
                     zero_filter=True, zero_filter_limit=0,
                     hundred_filter=True, hundred_filter_limit=0,
                     per_condition=False):
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
    per_condition: bool
        Must be True to filter quantity by condition individually. And False to filter by all samples.

    Returns
    -------
    df
        DataFrame without filtered positions
    """
    if df.iloc[:,:-2].empty:
        # print('The input is an empty DataFrame. In this case the function returns the input')
        return df


    if per_condition:
        array_filter = df.iloc[:,range(df.shape[1]-1)].columns # minus 1 to exclude diagnosis column
        for cond in df.iloc[:,-1].unique():
            # nan filter
            if nan_filter:
                array_filter = array_filter[ np.where( (df.loc[df.iloc[:,-1]==cond, array_filter ].isna().sum(axis=0))<=nan_filter_limit )[0]]
            # zero filter
            if zero_filter:
                array_filter = array_filter[ np.where( ((df.loc[df.iloc[:,-1]==cond, array_filter ]==0).sum(axis=0))<=zero_filter_limit )[0]]
            # hundred filter
            if hundred_filter:
                array_filter = array_filter[ np.where( ((df.loc[df.iloc[:,-1]==cond, array_filter ]==100).sum(axis=0))<=hundred_filter_limit )[0]]
        array_filter = array_filter.tolist()
        array_filter.append( df.columns[-1] )

    else:
        array_filter = df.iloc[:,range(df.shape[1]-2)].columns # minus 2 to exclude region column and diagnosis
        if nan_filter:
            # array_filter = array_filter[ np.where( (df.loc[:, array_filter ].isna().sum(axis=0))<=nan_filter_limit )[0]]
            aux = pd.DataFrame( df.loc[:, array_filter ].isna().sum(axis=0))
            array_filter = aux[aux[0]<=nan_filter_limit].index
        # zero filter
        if zero_filter:
            # array_filter = array_filter[ np.where( ((df.loc[:, array_filter ]==0).sum(axis=0))<=zero_filter_limit )[0]]
            aux = pd.DataFrame( (df.loc[:, array_filter ]==0).sum(axis=0)).sort_values(0)
            array_filter = aux[aux[0]<=zero_filter_limit].index
        # hundred filter
        if hundred_filter:
            # array_filter = array_filter[ np.where( ((df.loc[:, array_filter ]==100).sum(axis=0))<=hundred_filter_limit )[0]]
            aux = pd.DataFrame( (df.loc[:, array_filter ]==100).sum(axis=0)).sort_values(0)
            array_filter = aux[aux[0]<=zero_filter_limit].index

        array_filter = array_filter.tolist()
        array_filter.append( df.columns[-2] )
        array_filter.append( df.columns[-1] )
    return df[array_filter]

def replace_nan_by_zero(df, min_ratio_nan:float = 2/3):
    """
    Replace NaN values by zero  if the ratio of NaNs by total samples in one condition is greater than min_ratio_nan.

    Parameters
    ----------
    meta: df
        Frequency editing, Adenine or Guanine counts DataFrame of all samples analyzed. The DataFrame must be like merge_files output. The last two columns must be region and conditions.
    min_ratio_nan: float
        Must be a float representing the minimum ratio of NaN values of each condition that can be replaced by zero. Each column will be treated individually.

    Returns
    -------
    df
        DataFrame with NaNs values replaced by zeros
    """

    df_aux = df.copy()
    conditions, qnt_cond = np.unique(df_aux.iloc[:,-1].values, return_counts=True)
    cond_col = df_aux.columns[-1]
    for cond, qnt in zip(conditions, qnt_cond):
        aux = pd.DataFrame(df_aux[df_aux.diagnosis==cond].iloc[:,:-2].isna().sum()/qnt)
        cols = aux[aux[0]>=(min_ratio_nan)].index
        del(aux)
        for c in cols:
            df_aux.loc[df_aux[cond_col]==cond, c] = df_aux.loc[df_aux[cond_col]==cond, c].fillna(0)
    return df_aux

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
