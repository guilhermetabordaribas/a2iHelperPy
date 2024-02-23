import os
import pandas as pd
import anndata as ad
import numpy as np

def from_oligo(counts, obs=False, var=False, sep='\t', translation=False):
    adata = ad.AnnData( pd.read_csv(counts, sep=sep, index_col=0).T.apply(np.float32) )
    adata.var.index = adata.var.index.map(str)
    adata.obs.index = adata.obs.index.map(str)
    # print(adata)
    # print(adata.obs)
    # print(adata.var)
    if obs:
        # print(pd.read_csv(obs, sep=sep, index_col=0).head())
        obs = pd.read_csv(obs, sep=sep, index_col=0)
        obs.index = obs.index.map(str)
        adata.obs = adata.obs.merge( obs, right_index=True, left_index=True )
    if var:
        # print(pd.read_csv(var, sep=sep, index_col=0, header=None).head())
        var = pd.read_csv(var, sep=sep, index_col=0, header=None)
        var.index = var.index.map(str)
        adata.var = adata.var.merge( var, right_index=True, left_index=True )
    if translation:
        print('Translating', str(len(adata.var.index)), 'probes. It will take a long time. Consider pre-translate with the file https://github.guilhermetabordaribas....' )
        import mygene
        mg = mygene.MyGeneInfo()
        trans = [] # translated gene
        info = [] # information about translation process
        for r in adata.var.index:
            try:
                aux = mg.query( 'reporter:'+str(r) )
                if aux['hits']:
                    trans.append( str(aux['hits'][0]['symbol']) )
                    info.append('translation_ok')
                else:
                    trans.append( str(r) )
                    info.append("not_in_database")
            except:
                trans.append( str(r) )
                info.append("exception_was_found")

        adata.var['trans2Symbol'] = trans
        adata.var['trans_info'] = info
        adata.var.set_index('trans2Symbol', inplace=True)

    return adata.copy()

def from_HTseq_counts(path, obs=False, var=False, sep='\t', translation=False):
    count_list = os.listdir(path)
    df_list = []
    for file in count_list:
        df_list.append( pd.read_csv(os.path.join(path, file), sep=sep, header=None, names=['Genes',file]).set_index('Genes'))

    df_counts = pd.concat(df_list)
    print(df_counts)
    df_counts.set_index('Genes', inplace=True)

    adata = ad.AnnData( df_counts[~df_counts.index.str.startswith('__')].T.apply(np.float32).values )
    adata.obs = adata.obs.merge(df_counts[df_counts.index.str.startswith('__')].T, right_index=True, left_index=True )

    if obs:
        adata.obs = adata.obs.merge( pd.read_csv(obs, sep=sep, index_col=0), right_index=True, left_index=True )
    if var:
        adata.var = adata.var.merge( pd.read_csv(var, sep=sep, index_col=0), right_index=True, left_index=True )

    if translation:
        import mygene
        mg = mygene.MyGeneInfo()
        trans = [] # translated gene
        info = [] # information about translation process
        for r in adata.var.index:
            try:
                aux = mg.query( 'reporter:'+str(r) )
                if aux['hits']:
                    trans.append( str(aux['hits'][0]['symbol']) )
                    info.append('translation_ok')
                else:
                    trans.append( str(r) )
                    info.append("not_in_database")
            except:
                trans.append( str(r) )
                info.append("exception_was_found")

        adata.var['trans2Symbol'] = trans
        adata.var['trans_info'] = info
        adata.var.set_index('trans2Symbol', inplace=True)

    return adata.copy()
