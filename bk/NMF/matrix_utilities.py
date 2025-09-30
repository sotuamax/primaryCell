import numpy as np 
from scipy.stats import median_abs_deviation
from scipy import sparse 
from sklearn import preprocessing

def pair_pearson(X, Y):
    # feature correlation for X and Y
    xv = X - X.mean(axis = 0) # column mean
    yv = Y - Y.mean(axis = 0) # column mean 
    xvss = (xv*xv).sum(axis = 0) # column sum, variance
    yvss = (yv*yv).sum(axis = 0) # column sum, variance
    result = np.matmul(xv.transpose(), yv)/np.sqrt(np.outer(xvss, yvss))
    # np.clip(result, -1, 1)
    return result

def maxpcol(H):
    # scale H matrix to 0-1 range. this is based on column level, and the maximum value within each column are considered as the k-assignment for the m-sample in the kxm matrix
    # H_p = (H/H.sum(axis = 0))
    # H_score = 1+(1/np.log2(H.shape[0]))*np.multiply(H_p, np.log2(H_p+1))
    # np.nan_to_num(H_score, copy = False)
    # get the max for column
    # H_score = scale.fit_transform(H, axis = 0)
    sample_group = np.argmax(H, axis = 0)
    connectivity_m = np.array([i == j for i in sample_group for j in sample_group]).reshape((H.shape[1], H.shape[1])).astype(int)
    # max_value = np.max(H_scaled, axis = 0)
    scaler = preprocessing.MinMaxScaler()
    H_scaled = scaler.fit_transform(H)
    max_cut = np.median(H_scaled, axis = 0) + 2*median_abs_deviation(H_scaled, axis = 0)
    classfied = sum(max_cut<=1)/len(sample_group)
    return connectivity_m, classfied

def dispersion(consensus_m):
    # consensus matrix is the average connectivity matrix
    # return dispersion coefficient: =1 for perfect consensus matrix, [0,1) for scatter consensus matrix
    n = consensus_m.shape[1]
    dp = (1/n**2)*np.sum(4*(consensus_m - 0.5)**2)
    return dp

def maxprow(W, gene_index):
    scaler = preprocessing.MinMaxScaler()
    W_scaled = scaler.fit_transform(W.T).T
    # W_p = (W.T/W.sum(axis = 1)).transpose()
    # W_score = 1+1/np.log2(W.shape[1])*(W_p*np.log2(W_p))
    # np.nan_to_num(W_score, copy = False)
    gene_group = np.argmax(W_scaled, axis = 1)
    gene_downsample = gene_group[gene_index]
    #zeros = np.zeros((gene_group.size, gene_group.size))
    pos = np.array([(i,ii) for i,n in enumerate(gene_downsample) for ii,m in enumerate(gene_downsample) if n==m])#.reshape(gene_group.size, gene_group.size).astype("int64")
    # np.array([int(m==n) for n in gene_group for m in gene_group]).reshape(gene_group.size, gene_group.size)
    # np.put_along_axis(zeros, pos[:, 0].reshape(-1,), 1, pos[:,1].reshape(-1,))
    gene_m = sparse.csr_matrix((np.ones(len(pos)), (pos[:, 0].reshape(-1,), pos[:,1].reshape(-1,))), shape=(gene_downsample.size, gene_downsample.size))
    # max_value = np.max(W_score, axis = 1)
    feature_cut = np.median(W_scaled, axis = 1) + 3*median_abs_deviation(W_scaled, axis = 1)
    classfied = sum(feature_cut <= 1)/len(gene_group)
    # group_count = np.unique(W_maximum, return_counts = True)
    # to csr matrix 
    # gene_connect_sparse = sparse.csr_matrix(gene_connect)
    return gene_m, classfied

