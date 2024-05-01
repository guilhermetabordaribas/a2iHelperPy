import itertools
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from matplotlib.patches import Patch
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import seaborn as sns
import numpy as np
import pandas as pd
from adjustText import adjust_text
from statannotations.Annotator import Annotator
# from sklearn.decomposition import PCA
# from sklearn.manifold import TSNE

def boxplot(df, positions_to_plot:list = None, log_scale:bool = False,  ax=None, pvalue_list=None, figsize:tuple = None):
    if positions_to_plot == None:
        aux = df.drop(df.columns[-2], axis=1).melt(id_vars=df.columns[-1])
    else:
        if not isinstance(positions_to_plot, list):
            positions_to_plot = positions_to_plot.tolist()
        aux = df[positions_to_plot + [df.columns[-1]]].melt(id_vars=df.columns[-1])
    if log_scale:
        aux['value'] = np.log(aux.value+1)
    if ax == None:
        if figsize == None:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    x = 'variable'
    y = 'value'
    hue = df.columns[-1]

    hue_order = aux[aux.columns[0]].unique()
    order = aux.variable.unique()

    sns.boxplot(x=x, y=y, hue=hue, data=aux, ax=ax)
    ax.set_xlabel('Positions')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if log_scale:
        ax.set_ylabel('log(Editing Frequencies)')
    else:
        ax.set_ylabel('Editing Frequencies')
    ax.set_title(','.join(df.iloc[:,-2].unique()))

    if isinstance(pvalue_list, list):
        pairs = [tuple(itertools.product([p],hue_order)) for p in order]
        annot = Annotator(ax, pairs=pairs, x=x, y=y, hue=hue, order=order, hue_order=hue_order, data=aux)
        annot.configure(test=None, text_format='star', loc='inside', verbose=0)
        annot.set_pvalues(pvalue_list)
        ax, test_results = annot.annotate()

    return ax

def manhattanplot(df_or, df_pv, p_value_line:float = None, chr_order:list = [], ax=None, figsize:tuple = None):
    aux_m = pd.concat([df_or,df_pv]).iloc[:,:-2].T.reset_index()
    aux_m['-log10(pvalue)'] = -np.log10(aux_m['pvalue'])
    aux_m['index'] = aux_m['index'].astype(str)
    aux_m['chr'] = aux_m['index'].str.split('_').str[0]
    if ax == None:
        if figsize == None:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    x = 'chr'
    y = '-log10(pvalue)'
    hue = 'odds_ratio'
    order = None
    if chr_order:
        order = chr_order


    # sns.scatterplot(x=x, y=y, hue=hue, data=aux_m, ax=ax)
    sns.stripplot(x=x, y=y, hue=hue, order=order, palette='viridis_r', data=aux_m, ax=ax)
    ax.get_legend().remove()
    ax.set_xlabel('Chromosomes')
    ax.set_ylabel('-log10(pvalue)')
    ax.set_title(','.join(df_or.iloc[:,-1].unique())+' ('+','.join(df_or.iloc[:,-2].unique())+')')
    if p_value_line != None:
        ax.axhline(-np.log10(p_value_line), ls='--', lw=.5, color='black')

    cmap = mpl.cm.viridis_r
    norm = mpl.colors.Normalize(vmin=aux_m.odds_ratio.min(), vmax=aux_m.odds_ratio.max())

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=.15, pad=.01)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='Odds ratio', cax=cax)
    cbar.outline.set_color(None)

    return ax

def entropy_plot(entr, n_top:int = 50, log_scale:bool = False, ax=None, figsize:tuple = None ):
    aux = entr.melt(ignore_index=False).reset_index()
    aux = aux[ aux.variable.isin(aux[aux['index']==aux['index'].unique()[0]].sort_values('value', ascending=False).variable.values[:n_top]) ]
    order = aux[aux['index']==aux['index'].unique()[0]].sort_values('value', ascending=False).variable.values

    if log_scale:
        aux['value'] = np.log(aux.value)
    if ax == None:
        if figsize == None:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    x = 'variable'
    y = 'value'
    hue = 'index'

    sns.barplot(x=x, y=y, hue=hue,order=order, data=aux, ax=ax)
    ax.set_xlabel('Positions')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if log_scale:
        ax.set_ylabel('log(Entropy)')
    else:
        ax.set_ylabel('Entropy')
    ax.set_title(','.join(entr.iloc[:,-2].unique()))

    return ax

def corr_pearson_plot(p_corr, log_scale:bool = False,  p_value_line:float = None, ax = None, figsize:tuple = None):
    aux = p_corr.copy()
    x = 'Positions'
    y = 'pvalue'
    hue = 'R'
    if log_scale:
        aux['-log10(pvalue)'] = -np.log10(aux.pvalue)
        y = '-log10(pvalue)'
    if ax == None:
        if figsize == None:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    sns.scatterplot(x=x, y=y, hue=hue, data=aux, ax=ax)
    ax.get_legend().remove()
    ax.set_xlabel('Positions')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    if log_scale:
        ax.set_ylabel('-log10(pvalue)')
        if p_value_line != None:
            ax.axhline(-np.log10(p_value_line), ls='--', lw=.5, color='black')
    else:
        ax.set_ylabel('pvalue')
        if p_value_line != None:
            ax.axhline(p_value_line, ls='--', lw=.5, color='black')
    ax.set_title(','.join(df.iloc[:,-2].unique()))

    cmap = mpl.cm.viridis_r
    norm = mpl.colors.Normalize(vmin=aux.R.min(), vmax=aux.R.max())
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size=.15, pad=.01)
    cbar = plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, label='R', cax=cax)
    cbar.outline.set_color(None)

    return ax

def volcanoplot(data, pv_col='padj', pv_lim=0.1, logFC_col='log2FoldChange', logFC_lim=1.5, gene_col=False, figsize=None, ax=None, use_adjusttext=False, text=None):
    df = data.copy()
    df.dropna(inplace=True)
    df['-log10(padj)'] = -np.log10(df[pv_col])
    df['log2FoldChange_abs'] = df[logFC_col].abs()

    condlist = [(df[logFC_col] >= logFC_lim) & (df[pv_col] <= pv_lim), (df[logFC_col] <= -logFC_lim) & (df[pv_col] <= pv_lim)]
    choicelist = ['Upregulated','Downregulated']

    df['regulation'] = np.select(condlist, choicelist, default='Not-significant')

    if ax == None:
        if figsize:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    hue = 'regulation'
    hue_order = ['Upregulated', 'Downregulated']#,'Not-significant']
    palette = ['red', 'royalblue']#, 'grey']
    x = logFC_col
    y = '-log10(padj)'

    # sns.scatterplot(x=x, y=y, hue=hue, hue_order=hue_order, palette=palette, edgecolor=None, data=df, ax=ax)
    sns.scatterplot(x=x, y=y, color='gray', edgecolor=None, alpha=.25, data=df[df['regulation']=='Not-significant'], ax=ax)
    # print('hey')
    sns.scatterplot(x=x, y=y, hue=hue, hue_order=hue_order, palette=palette, edgecolor='black', data=df[df['regulation']!='Not-significant'], ax=ax)
    ax.axvline(logFC_lim, color='black', ls='dashed', lw=1.5, alpha=0.6)
    ax.axvline(-logFC_lim, color='black', ls='dashed', lw=1.5, alpha=0.6)
    ax.axhline(-np.log10(pv_lim), color='black', ls='dashed', lw=1.5, alpha=0.6)

    if use_adjusttext:
        if text == None:
            texts = []
            if gene_col:
                list_ = df[df.regulation!='Not-significant'].sort_values(['-log10(padj)','log2FoldChange_abs'], ascending=False)[[logFC_col, '-log10(padj)', gene_col]].values[:30]
            else:
                list_ = df[df.regulation!='Not-significant'].reset_index().sort_values(['-log10(padj)','log2FoldChange_abs'], ascending=False)[[logFC_col, '-log10(padj)','index']].values[:30]
            for x, y, l in list_:
                if False:# l in ['CTSS','TAP2','TAP1','HLA-DPA1']:
                    texts.append(ax.text(x, y, l, fontsize=lgnd_fz, bbox={'facecolor':'white', 'alpha':.5, 'edgecolor':'white', 'pad':0}))#
                else:
                    texts.append(ax.text(x, y, l, fontsize=8))#, bbox={'facecolor':'white', 'alpha':.5, 'edgecolor':'white'}

        else:
            texts = []
            if gene_col:
                list_ = df.loc[text, :][[logFC_col, '-log10(padj)', gene_col]]
            else:
                list_ = df.loc[text, :][[logFC_col, '-log10(padj)','index']]
            for x, y, l in list_:
                texts.append(ax.text(x, y, l, fontsize=8))
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color='k', lw=1.1),
                    expand_text=(1.2, 1.2), expand_points=(1.2, 1.2),
                    precision=0.0001, lim=5000
                   )

    return ax

def pca(data, condition, hue=False, ax=None, figsize=None, conf_ellipse=False):
    X = data.X

    pca = PCA(n_components=2)
    X = pca.fit_transform(X)
    df_pca = pd.DataFrame(X)

    if ax == None:
        if figsize:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    if hue:
        hue='Type'
    else:
        hue = None

    m = data.obs[condition].values
    df_pca['Type'] = m
    sns.scatterplot(x=0, y=1, hue=hue, data=df_pca, ax=ax)

    if conf_ellipse:
        for tl in df_pca['Type'].unique():
            confidence_ellipse(df_pca[df_pca.Type==tl][0].values,
                               df_pca[df_pca.Type==tl][1].values, ax, edgecolor='#000000', alpha=.7, n_std=3)

    var = pca.explained_variance_ratio_

    ax.set_xlabel('PC1 ('+str(round(var[0]*100, 2))+')%' )
    ax.set_ylabel('PC2 ('+str(round(var[1]*100, 2))+')%' )

def tsne(data, condition, hue=False, ax=None, figsize=None, conf_ellipse=False):
    X = data.X

    tsne = TSNE(n_components=2)
    X = tsne.fit_transform(X)
    df_pca = pd.DataFrame(X)

    if ax == None:
        if figsize:
            f, ax = plt.subplots()
        else:
            f, ax = plt.subplots(figsize=figsize)

    if hue:
        hue='Type'
    else:
        hue = None

    m = data.obs[condition].values
    df_pca['Type'] = m
    sns.scatterplot(x=0, y=1, hue=hue, data=df_pca, ax=ax)

    if conf_ellipse:
        for tl in df_pca['Type'].unique():
            confidence_ellipse(df_pca[df_pca.Type==tl][0].values,
                               df_pca[df_pca.Type==tl][1].values, ax, edgecolor='#000000', alpha=.7, n_std=3)


def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`

    Returns
    -------
    matplotlib.patches.Ellipse
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensional dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2,
                      facecolor=facecolor, **kwargs)

    # Calculating the standard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the standard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
