def transform_data(count_df):
    """provide data in n_samples x n_features, and n_samples are of interest for clustering analysis"""
    dhs_log = np.log2(count_df+1) # to alleviate the effects from outlier
    return dhs_log  

def collect_cluster(labels, block_size, size, output):
    n_clusters = max(labels)+1
    start = 0
    end = 0
    W_bins = {i:[] for i in range(n_clusters)}
    for i in range(size):
        # Get labels for each partial matrix
        end += block_size[i]
        proc_labels = labels[start:end]
        start = end
        # Add parts of matrix to full table
        W_part = np.load(f"{output}_W_{i}.npy")
        # W_part.columns = range(W_part.shape[1])
        for i,lab in enumerate(proc_labels):
            if lab != -1:
                W_bins[lab].append(W_part[:,i])
    return W_bins

def cluster_centroid(W_cluster, count_df):
    W_final = pd.DataFrame(columns=range(len(W_cluster)),index=count_df.index)
    for lab in range(len(W_cluster)):
        W_clust = W_cluster[lab]
        # Add in rest of components
        # Add centroid of cluster to final S matrix
        W_final[lab] = np.array(W_clust).T.mean(axis=1)
    return W_final