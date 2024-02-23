from sklearn.base import BaseEstimator, TransformerMixin
import numpy as np

class TPM(BaseEstimator, TransformerMixin):
    def __init__(self, gene_len, idx='all', scaler=False):
        self.gene_len = gene_len
        self.idx = idx
        self.scaler = scaler
        self.not_zero_sum = ''
        # print('\n>>>init() called.\n')
    def fit(self, X, y=None):
        # print('\n>>>fit() called.\n')
        return self
    def transform(self, X, y=None):
        # print('\n>>>transform() called.\n')
        X_ = X.copy()
        a = X_ / (self.gene_len) / 1e3
        # a = a[:, np.where(a.sum(axis=0)>10)[0] ]
        np.seterr(all='warn')
        tpm = ( a.T / (a.sum(axis=1) / 1e6) ).T #sum of each replicate
        # tpm = a.T * 1e6
        if self.scaler == 'log':
            tpm = np.log2(tpm+1)
        elif self.scaler != False:
            tpm = self.scaler.fit_transform(tpm)
        if  isinstance(self.idx, str):
            if self.idx == 'all':
                return tpm
        else:
            return tpm[:, self.idx]

class SCALE_FACTOR(BaseEstimator, TransformerMixin):
    def __init__(self, scaler=False):
        self.scaler = scaler
        self.not_zero_sum = ''
        # print('\n>>>init() called.\n')
    def fit(self, X, y=None):
        # print('\n>>>fit() called.\n')
        return self
    def deseq2_scl_fac(self, X):
        # t = X.copy()
        t = np.log(X.copy())

        m = np.mean(t, axis=0)
        t = t[:,np.where(np.isfinite(m))[0]]
        m = m[np.isfinite(m)]

        median = np.median(t-m, axis=1)

        if self.scaler != False:
            return self.scaler.fit_transform( (X.copy().T/np.e**median).T )
        else:
            return (X.copy().T/np.e**median).T
