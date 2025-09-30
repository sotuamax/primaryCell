#!/usr/bin/env python3
from mpi4py import MPI
import argparse 
import logging 
import pandas as pd 
import numpy as np 
from sklearn.metrics.pairwise import pairwise_distances
from sklearn.cluster import DBSCAN
from sklearn.decomposition import NMF 
from sklearn.metrics import silhouette_score
from scipy import sparse 
import os 
from NMF.matrix_utilities import *
from utility import timeit 
import time
import sys
from itertools import chain 

def args_parser():
    '''parser the argument from terminal command'''
    parser=argparse.ArgumentParser(prog = "PROG", formatter_class = argparse.RawDescriptionHelpFormatter, description="")
    parser.add_argument("count", help = "count file")
    parser.add_argument("-k", "--components", help = "number of components", type = int, required = True)
    parser.add_argument("-r", "--repeats", help = "number of repeat runs ", type = int, default = 1)
    parser.add_argument("-app", "--appendt", default = 0, type = int, help = "naming ritual append of last run")
    parser.add_argument("-transform", "--transform", action = "store_true", help = "transform data to log2 scale")
    parser.add_argument("-gene_list", "--gene_list", required = False, help = "down sample gene list")
    parser.add_argument("-consensus", "--consensus", default = 0, type = int, help = "report consensus for multiple repeat runs (integer)")
    parser.add_argument("-dist", "--dist", action = "store_true", help = "distance for independent runs.")
    parser.add_argument("-skip", "--skip", action = "store_true", help = "skip decomposition step")
    parser.add_argument("-dir", "--dir", help = "output directory")
    parser.add_argument("-o", "--output", help = "output name")
    args=parser.parse_args()
    return args

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter(action = "ignore")
    args = args_parser()
    count = args.count
    output = os.path.join(args.dir, args.output)
    repeats = args.repeats
    components = args.components
    count_df = pd.read_table(count, sep = "\t", header = 0, index_col = 0)
    # create directory for outputs to store
    if rank == 0:
        start = time.time()
        if not os.path.exists(args.dir):
            os.mkdir(args.dir)
        test = 1
    else:
        test = 0
    test = comm.bcast(test, root = 0)
    # append logging info for file
    logging.basicConfig(filename = output+".log", filemode = "a", format = '%(asctime)s  %(message)s', datefmt = "%H:%M:%S", level = logging.INFO)
    # not skip, then perform repeat runs 
    if not args.skip:
        if args.transform:
            if rank == 0:
                logging.info("Tranform count data")
            count_df = transform_data(count_df)
            count_df = count_df[(count_df < 1).sum(axis = 1) != count_df.shape[1]]
            count_df[count_df < 1] = 0
            if rank == 0: 
                pd.DataFrame(count_df.index).to_csv(f"{output}_index.txt", sep = "\t", header = False, index = False)
        if rank == 0:
            worker_tasks = {w:[] for w in range(size)}
            w_idx = 0
            for r in range(args.appendt, repeats+args.appendt):
                worker_tasks[w_idx].append(r)
                w_idx = (w_idx + 1) % size
        else:
            worker_tasks = None 
        worker_tasks = comm.bcast(worker_tasks, root = 0)
        ### repeat run 
        if rank == 0:
            logging.info(f"Start NMF with {components} components")
            logging.info(f"Repeat from {args.appendt} - {repeats+args.appendt}")
        for r in worker_tasks[rank]:
            nmf = NMF(n_components = components)
            W = nmf.fit_transform(count_df)
            H = nmf.components_
            err = nmf.reconstruction_err_
            np.savez(f"{output}_{r}.npz", W=W, H=H, err=np.array([err]))
            # sparse.save_npz(f"{output}_{r}_sparse.npz", clust)
            logging.info(f"Finished repeat {r}")
        del worker_tasks
    # when call consensus, then perform consensus call from repeat runs
    if args.consensus != 0:
        if rank == 0:
            logging.info(f"Consensus call for range 0 to {args.consensus}")
            worker_tasks = {w:[] for w in range(size)}
            w_idx = 0
            for r in range(args.consensus):
                worker_tasks[w_idx].append(r)
                w_idx = (w_idx + 1) % size
        else:
            worker_tasks = None
        worker_tasks = comm.bcast(worker_tasks, root = 0)
        err_list = list(); hclust_list = list(); wclust_list = list()
        gene_index = None
        if args.gene_list != None:
            gene_list = pd.read_table(args.gene_list, sep = "\t", header = None, index_col = 0)
            gene_index = np.argwhere(count_df.index.isin(gene_list.index))
        for r in worker_tasks[rank]:
            repeat_handle = np.load(f"{output}_{r}.npz")
            hclust,h_classifed = maxpcol(repeat_handle["H"])
            if not gene_index is None:
                wclust,w_classified = maxprow(repeat_handle["W"], gene_index)
                wclust_list.append(wclust.toarray())
                logging.info(f"classified W {w_classified}; H {h_classifed}")
            logging.info(f"read repeat {r}")
            logging.info(f"classified H {h_classifed}")
            err_list.append(repeat_handle["err"])
            hclust_list.append(hclust)
        hclust_sum = sum(hclust_list)
        del hclust_list; 
        hclust_all = comm.gather(hclust_sum, root = 0)
        del hclust_sum; #del wclust_sum
        if not gene_index is None:
            wclust_sum = sum(wclust_list)
            wclust_all = comm.gather(wclust_sum, root = 0)
            del wclust_list
        err_all = comm.gather(err_list, root = 0)
        if rank == 0:
            logging.info("Integrate repeats for divergence calculation")
            h_dispersion = dispersion(sum(hclust_all)/args.consensus)
            # wclust_all = list(chain.from_iterable(wclust_all))
            err_all = list(chain.from_iterable(err_all))
            # err_np.mean(err_list)
            err_mean = round(np.mean(err_all), 2)
            err_df = pd.DataFrame([h_dispersion, err_mean], index = ["h_dispersion", "error"])
            if not gene_index is None:
                w_dispersion = dispersion(sum(wclust_all)/args.consensus)
                logging.info(f"Dispersion value on gene group {w_dispersion}")
                w_err = pd.DataFrame([w_dispersion], index = ["w_dispersion"])
                err_df = pd.concat([w_err, err_df], axis = 0)
            err_df.to_csv(f"{output}_err.txt", sep = "\t", index = True, header = False)
            logging.info(f"reconstrution error (mean) {err_mean}")
            logging.info(f"Dispersion value on sample group {h_dispersion}")
    if args.dist != None:
        dist_pair = [(n,m) for n in range(repeats) for m in range(n, repeats)]
        
        dist_dict = dict();dist_all = list()
        for n in range(repeats):
            dist_block = list()
            for m in range(repeats):
                if (n,m) in dist_pair:
                    d1 = np.load(f"{output}_{n}.npz")["W"]
                    d2 = np.load(f"{output}_{m}.npz")["W"]
                    dist_dict[(n,m)] = pair_pearson(d1, d2)
                    dist_block.append(dist_dict[(n,m)])
                else:
                    dist_block.append(dist_dict[(m,n)].T)
            dist_all.append(dist_block)
        dist_mat = np.block(dist_all); dist_mat = np.clip(dist_mat, -1, 1)
        dbscan = DBSCAN(eps=0.2, min_samples=5, metric='precomputed') # label as dhs_log.index
        labels = dbscan.fit_predict(1-abs(dist_mat)) 
        lab_dict = {r:[] for r in range(labels.max()+1)}
        for index, lab in enumerate(labels):
            if lab != -1:
                lab_dict[lab].append(index)


    if rank == 0:
        logging.info(timeit(start))
        logging.info("Finished!")
    exit(0)

if __name__ == "__main__":
    main()


