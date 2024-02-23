import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from matplotlib.patches import Patch
import seaborn as sns
import numpy as np
import pandas as pd
# from adjustText import adjust_text
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE


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
