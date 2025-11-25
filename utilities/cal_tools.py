import numpy as np 
from scipy.stats import rankdata
import pandas as pd 

def entropy(mat, static = False):
    """
    from https://www.nature.com/articles/s41586-020-2151-x
    interpretation: 
    higher value indicates a variable profile.
    lower value indicates a uniform profile; 
    """
    mat_index = mat.index
    mat = np.array(mat)
    # get the relative gene expression for each sample 
    frac_mat = np.array(mat)/(mat.sum(axis = 1).reshape(-1,1))#[:, None]
    # null model: each sample has equal expression level at per gene level
    bg_mat = np.ones(mat.shape)
    null_mat = bg_mat/(bg_mat.sum(axis = 1).reshape(-1,1))#[:, None]
    # calculate gene entropy
    rel_entropy = (frac_mat*np.log2(frac_mat/null_mat+1)).sum(axis = 1)
    assert rel_entropy.shape[0] == mat.shape[0], print("return results is not match to input dimension")
    if static:
        rel_entropy = 1/rel_entropy * mat.mean(axis = 1)
    rel_entropy = pd.DataFrame({"entropy":rel_entropy}, index = mat_index)
    return rel_entropy

def CoeVar(mat):
    """
    also known as relative standard deviation (RSD) - a standardized measure of the dispersion of a probability distribution. 
    the extent of variability in relation to the mean of the population. 
    - CV=σ/μ
        where: σ=standard deviation; μ=mean
    interpretation: 
    higher value indicates a variable profile.
    lower value indicates a uniform profile; 
    """
    cov_df = pd.DataFrame({"CoeVar":mat.std(axis = 1)/mat.mean(axis = 1)}, index = mat.index)
    return cov_df

def rowwise_spearman(X, Y):
    """
    Row-wise Spearman correlation using average ranks (tie-aware)
    X.shape - (a, b), Y.shape - (a, b)
    Taking each row of X and Y in a paired fashion, calculate Spearman's correlation for each row. 
    return shape - (a, 1)
    """
    # Rank each row with ties handled as average
    Xr = np.apply_along_axis(rankdata, 1, X, method='average')
    Yr = np.apply_along_axis(rankdata, 1, Y, method='average')
    # Compute row-wise Pearson on ranks
    Xc = Xr - Xr.mean(axis=1, keepdims=True)
    Yc = Yr - Yr.mean(axis=1, keepdims=True)
    num = np.sum(Xc * Yc, axis=1)
    den = np.sqrt(np.sum(Xc**2, axis=1) * np.sum(Yc**2, axis=1))
    return num / den

def pair_pearson(X, Y):
    """
    from https://cancerdatascience.org/blog/posts/pearson-correlation/ 
    X.shape - X.shape - (a, b), Y.shape - (a, b)
    X.shape - (a, b), Y.shape - (a, b)
    Taking each row of X and Y in all combination possibilities, calculate Pearson's correlation 
    return shape - (a, a)
    """
    xv = X - X.mean(axis = 0) # column mean
    yv = Y - Y.mean(axis = 0) # column mean 
    xvss = (xv*xv).sum(axis = 0) # column sum, variance
    yvss = (yv*yv).sum(axis = 0) # column sum, variance
    result = np.matmul(xv.transpose(), yv)/np.sqrt(np.outer(xvss, yvss))
    return result

