import numpy as np 
from scipy.stats import rankdata
import pandas as pd 
from scipy.stats import chi2_contingency, fisher_exact

def quantile_normalize_np(x):
    from sklearn.preprocessing import quantile_transform
    """
    Quantile normalize a 2D numpy array (rows = features, cols = samples).
    """
    x = np.asarray(x)
    ranks = np.argsort(np.argsort(x))
    return ranks / (len(x) - 1)


def local_linear_regression(array, window=2):
    """
    Docstring for local_linear_regression: calculate the local slope for 1-D array

    :param array: array used to calculate local slopes
    :param window: window size to use for slope
    :param n: number of threads to use 
    """
    slopes = np.full(len(array), np.nan)
    n = len(array)
    x = np.arange(-window, window+1).reshape(-1,1)
    denom = (x.T @ x)[0, 0]
    for i in range(window, n - window):
        y = array[i - window : i + window + 1]
        slopes[i] = (x[:, 0] @ y) / denom
    return slopes

def derivative(arr, side = 3):
    """
    Caculate derivative of array by taking the side values for their difference (centered on one value). 
    """
    arr = np.asarray(arr, dtype=float)
    kernel = np.concatenate([
        np.ones(side) / side,
        np.zeros(1),
        -np.ones(side) / side
    ]) # convolve flips the kernel internelly
    return np.convolve(arr, kernel, mode="same") 

def distance_from_median(df):
    """
    Distance (absolute deviations) from the median per row
    """
    row_median = np.median(df, axis = 1)
    diff_median = df.sub(row_median, axis = 0).abs()
    return diff_median

def distance_from_mean(df):
    """
    distance from the mean per row
    """
    row_mean = df.mean(axis=1)
    # element-wise distance from row mean
    distance = df.sub(row_mean, axis=0).abs()
    return distance

def fisher_enrich(df, padjust=True, as_table=False, return_dict = True):
    """
    df: DataFrame with two categorical columns OR
        a contingency table (if as_table=True)

    Returns dictionary with:
      observed, expected, residual, odds.ratio,
      p, CI.lower, CI.upper, padj (optional)
    """
    # ---- build contingency table ----
    if not as_table:
        df = df.copy()
        row_feat = df.iloc[:, 0].astype("category")
        col_feat = df.iloc[:, 1].astype("category")
        tbl = pd.crosstab(row_feat, col_feat)
    else:
        tbl = df.copy()
    # ---- chi-square independence test ----
    chi2, pval, dof, expected = chi2_contingency(tbl)
    if pval >= 0.01:
        print("Two variables are independent (H0).")
        return None
    else:
        print("Two variables are dependent.")
        print("Perform fisher's exact test ...")
    observed = tbl.values
    expected = expected
    residual = (observed - expected) / np.sqrt(expected)
    rowCats = tbl.index
    colCats = tbl.columns
    # ---- containers ----
    pTbl = pd.DataFrame(index=rowCats, columns=colCats, dtype=float)
    orTbl = pd.DataFrame(index=rowCats, columns=colCats, dtype=float)
    lowerci = pd.DataFrame(index=rowCats, columns=colCats, dtype=float)
    upperci = pd.DataFrame(index=rowCats, columns=colCats, dtype=float)
    # ---- Fisher's exact test for each rc, cc combination ----
    for rc in rowCats:
        for cc in colCats:
            rclass = (df.iloc[:, 0] == rc)
            cclass = (df.iloc[:, 1] == cc)
            # 2x2 contingency table
            table_2x2 = np.array([
                [(rclass & cclass).sum(), (rclass & ~cclass).sum()],
                [(~rclass & cclass).sum(), (~rclass & ~cclass).sum()]
            ])
            # Fisher exact test
            odds_ratio, p = fisher_exact(table_2x2, alternative="two-sided")
            # CI (scipy does not provide CI; compute via Woolf log method)
            # Avoid log issues
            a, b, c, d = table_2x2.flatten()
            if 0 in [a, b, c, d]:
                # Haldane-Anscombe correction
                a += 0.5; b += 0.5; c += 0.5; d += 0.5
            # log(OR) CI
            se = np.sqrt(1/a + 1/b + 1/c + 1/d)
            log_or = np.log((a*d) / (b*c))
            ci_low = np.exp(log_or - 1.96 * se)
            ci_up = np.exp(log_or + 1.96 * se)
            pTbl.loc[rc, cc] = p
            orTbl.loc[rc, cc] = odds_ratio
            lowerci.loc[rc, cc] = ci_low
            upperci.loc[rc, cc] = ci_up
    # ---- padjust ----
    enrich_dict = {
        "observed": pd.DataFrame(observed, index=rowCats, columns=colCats),
        "expect": pd.DataFrame(expected, index=rowCats, columns=colCats),
        "residual": pd.DataFrame(residual, index=rowCats, columns=colCats),
        "odds.ratio": orTbl,
        "p": pTbl,
        "CI.lower": lowerci,
        "CI.upper": upperci,
    }
    if padjust:
        from statsmodels.stats.multitest import multipletests
        flat_p = pTbl.values.flatten()
        padj = multipletests(flat_p, method="fdr_bh")[1]
        padjTbl = pd.DataFrame(
            padj.reshape(pTbl.shape),
            index=rowCats, columns=colCats
        )
        enrich_dict["padj"] = padjTbl
    if return_dict:
        return enrich_dict 
    else:
        enrich_df = enrich_dict["odds.ratio"].reset_index().melt(id_vars=enrich_dict["odds.ratio"].index.name, var_name="feature", value_name="oddsRatio")
        enrich_p = enrich_dict["p"].reset_index().melt(id_vars=enrich_dict["p"].index.name, var_name="feature", value_name="p")
        enrich_df = pd.merge(enrich_df, enrich_p, on = enrich_df.columns[:2].tolist())
        if padjust:
            enrich_padj = enrich_dict["padj"].reset_index().melt(id_vars=enrich_dict["padj"].index.name, var_name="feature", value_name="padj")
            enrich_df = pd.merge(enrich_df, enrich_padj, on = enrich_df.columns[:2].tolist())
        return enrich_df


def pipe_specific(df, group_df, select_group):
    """
    Docstring for pipe_specific
    
    :param df: data frame in a matrix format 
    :param group_df: group dataframe, each sample match to a column in df
    :param select_group: select specific group for specific test (t-test, z-test)
    """
    assert sorted(group_df.index) == sorted(df.columns), print("group is not match to data.")
    group_df = group_df.reindex(df.columns, copy = True)
    assert group_df.index.tolist() == df.columns.tolist(), print("group is not properly ordered in accordance with data.")
    group_df.columns = ["group"]
    assert select_group in group_df["group"].tolist(), print("assigned group does not present in groups.")
    if len(group_df.query("group == @select_group")) > 1:
        group_stat = t_statistic(df, group_df["group"], select_group) 
    else: 
        group_stat = z_statistic(df, group_df["group"], select_group)
    return group_stat 

def t_statistic(df:pd.DataFrame, group:np.array, select_group:str):
    """
    vectorized t-tests to accelerate on each row of df, by comparing select_group to others. 

    # apply to each row (vectorized version)
    x <- v[grp == "A"]; y <- v[grp == "B"]
    tt <- t.test(x, y, alternative = "greater")
    #m <- length(x); n <- length(y)
    # t <- tt$statistic[[1]] # rank sum for group A
    # tt$stderr for how variability in the difference between x and y 
    c(statistic = tt$statistic[[1]], p = tt$p.value, contactA = mean(x), contactB = mean(y), standard_error= tt$stderr)
    """ 
    from scipy.stats import t, norm
    mat = np.array(df)
    assert mat.shape[1] == len(group), print("Given data match to group in length")
    groupA = group == select_group
    groupB = group != select_group
    A = mat[:, groupA]; B = mat[:, groupB]
    nA = A.shape[1]; nB = B.shape[1]
    meanA = A.mean(axis = 1); meanB = B.mean(axis = 1); diff = meanA - meanB
    varA = A.var(axis = 1); varB = B.var(axis = 1) # by default, ddof=1, N - ddof is delta degree of freedom 
    stderr = np.sqrt(varA / nA + varB / nB) # this value differs from default R 
    t_stat = diff / stderr
    t_stat[np.isnan(t_stat)] = 0
    df_num = (varA / nA + varB / nB) ** 2
    df_den = (varA**2) / (nA**2 * (nA - 1)) + (varB**2) / (nB**2 * (nB - 1))
    df_eff = df_num / df_den
    #sd_pooled = np.sqrt(((nA - 1) * varA + (nB - 1) * varB) / (nA + nB - 2))
    #effect_size = diff / sd_pooled
    #alpha = 0.05
    #t_crit = t.ppf(1 - alpha/2, df_eff)   # two-sided critical value
    #CI_low  = diff - t_crit * stderr
    #CI_high = diff + t_crit * stderr
    pvals = 1 - t.cdf(t_stat, df_eff)
    if nA == 1:
        pvals = 1 - norm.cdf(t_stat)
    out = pd.DataFrame({
        "statistic": t_stat,
        "A": meanA,
        "B": meanB,
        "SE": stderr, # standard error 
        "df":df_eff,
        #"CI_low": CI_low, 
        #"CI_high": CI_high,
        #"effect": effect_size, 
        "p": pvals # pvalue 
    }, index=df.index)
    return out

def z_statistic(df:pd.DataFrame, group:np.array, select_group:str):
    """
    when group has one sample only. 
    R version: 
    ztest_row <- function(v, grp = groups) {
    x <- as.numeric(v[grp == "A"]); y <- v[grp == "B"]
    mu <- mean(y); sigma <- sd(y)
    z <- (x - mu)/sigma 
    p_value <- 1 - pnorm(z)
    # tt$stderr for how variability in the difference between x and y 
    c(statistic = z, p = p_value, contactA = mean(x), contactB = mu, error= sqrt(sigma))
    }
    """
    from scipy.stats import norm
    mat = np.array(df)
    assert mat.shape[1] == len(group), print("Given data match to group in length")
    groupA = group == select_group
    groupB = group != select_group
    assert sum(groupA) == 1, print("selected group has more than 1 sample, use t-statistic instead.")
    A = mat[:, groupA]; B = mat[:, groupB]
    nA = A.shape[1]; nB = B.shape[1]
    meanA = A.mean(axis = 1); meanB = B.mean(axis = 1)
    varA = A.var(axis = 1); varB = B.var(axis = 1)
    stderr = np.sqrt(varA / nA + varB / nB)
    z_stat = (meanA-meanB)/stderr
    z_stat[np.isnan(z_stat)] = 0
    pvals = 1 - norm.cdf(z_stat) # cdf is cumulative distribution function that a random variable X is less than x
    out = pd.DataFrame({
        "statistic": z_stat,
        "A": meanA,
        "B": meanB,
        "SE": stderr, 
        "p": pvals,
        
    }, index=df.index)
    return out

def entropy(mat, static = False):
    """
    from https://www.nature.com/articles/s41586-020-2151-x
    the differences between observed fraction / expected fraction of values in all variable values
    This value is not propriate for matrix with negative values. 
    
    interpretation: 
    higher value indicates a variable profile.
    lower value indicates a conserved profile; 

    Return: 
    entropy in a dataframe 
    """
    mat_index = mat.index
    mat = np.array(mat)
    # get the relative gene expression for each sample (observed fraction)
    frac_mat = np.array(mat)/(mat.sum(axis = 1).reshape(-1,1))#[:, None]
    # null model: each sample has equal expression level at per gene level
    bg_mat = np.ones(mat.shape)
    null_mat = bg_mat/(bg_mat.sum(axis = 1).reshape(-1,1))#[:, None] (expected fraction)
    # calculate gene entropy (differences in observed fraction / expected fraction)
    rel_entropy = (frac_mat*np.log2(frac_mat/null_mat+1)).sum(axis = 1) # note: log2(+1) was added by me, to avoid 0 in the calculation
    assert rel_entropy.shape[0] == mat.shape[0], print("return results is not match to input dimension")
    rel_entropy_df = pd.DataFrame({"entropy":rel_entropy}, index = mat_index)
    if static:
        static_value = 1/rel_entropy * mat.mean(axis = 1)
        rel_entropy_df["static"] = static_value
    return rel_entropy_df

def CoeVar(mat):
    """
    It is not appropriate for data with both negative and positive values (mean may be close to 0, and normalized by mean can be over-valued)
    also known as relative standard deviation (RSD) - a standardized measure of the dispersion of a probability distribution. 
    the extent of variability in relation to the mean of the population. 
    - CV=σ/μ
        where: σ=standard deviation; μ=mean
    interpretation: 
    higher value indicates a variable profile.
    lower value indicates a uniform profile; 
    return: 
    Coefficient Variance in dataframe 
    """
    cov_df = pd.DataFrame({"CoeVar":mat.std(axis = 1)/mat.mean(axis = 1)}, index = mat.index)
    return cov_df

def row_min(mat):
    min_df = pd.DataFrame({"Min":mat.min(axis = 1)}, index = mat.index)
    return min_df 

def row_max(mat):
    max_df = pd.DataFrame({"Max":mat.max(axis = 1)}, index = mat.index)
    return max_df 

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

