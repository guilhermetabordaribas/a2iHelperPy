import numpy as np
# from boruta import BorutaPy #for boruta, there is a error with numpy.int. I need to replace np.int by np.int64 or np.int32. the samething to np.float and np.bool to np.bool_
import numbers

def boruta_selection(data, model, class_col='target', scaler=None, random_state=None):
    X = data.X # insert option to layers
    y = data.obs[class_col]

    if scaler != None:
        X_scaler = scaler.fit_transform(X)
    else:
        X_scaler = X

    feat_selector = BorutaPy(model, n_estimators='auto', verbose=0, random_state=random_state)
    feat_selector.fit(X_scaler, y)

    feature_ranks = list(zip(data.var.index,
                         feat_selector.ranking_,
                         feat_selector.support_))

    # iterate through and print out the results

    result = {'genes':[], 'rank':[], 'support':[]}
    for feat in feature_ranks:
        if feat[2]:
            result['genes'].append(feat[0])
            result['rank'].append(feat[1])
            result['support'].append(feat[2])
            # print('Feature: {:<25} Rank: {},  Keep: {}'.format(feat[0], feat[1], feat[2]))
    return result

def set_boruta_selection(data, model, class_col='target', scaler=None, random_state=None, n_set=5, sample_size=False, replace=False):
    results = []
    for i,set in enumerate(set_balance_resample(data.obs[class_col].values, n_set=n_set, random_state=random_state, sample_size=sample_size, replace=replace)):
        results.append(boruta_selection(data, model, class_col=class_col, scaler=scaler, random_state=random_state+i))

    return results

def balance_resample(y_var, random_state=42, sample_size=False, replace=False):
    v, c = np.unique(y_var, return_counts=True)
    if sample_size == False:
        sample_size = c.min()
    else:
        if sample_size > c.min():
            print("Warning, sample_size greather than the minimum class counts. Truning sample_size =", str(c.min()))
            sample_size = c.min()

    idx = np.array([], dtype=int)
    for v_ in v:
        np.random.seed(random_state)
        idx = np.concatenate( ( idx, np.random.choice(np.where(y_var == v_)[0], size=sample_size, replace=replace) ))

    return idx

def set_balance_resample(y_var, n_set=5, random_state=42, sample_size=False, replace=False):
    set_idx = []
    for n in range(n_set):
        set_idx.append( balance_resample(y_var, random_state+n, sample_size, replace) )

    return set_idx

# def
