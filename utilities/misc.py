#!/usr/bin/env python3
import pandas as pd
import numpy as np 
import pysam 
import time 
import os 
import subprocess 


def ignore_warning():
    """
    Ignore warnings 
    """
    import sys 
    if not sys.warnoptions:
        import warnings
        warnings.simplefilter("ignore")

def module_load(modules:list):
    """
    Load modules and set up new environment path
    Input: 
    modules: list of modules to load 
    Output:
    new environment path (if using subprocess, set up env as the new environment to be able to locate the path)
    """
    env = os.environ.copy()
    module = " ".join(modules)
    # print(module)
    
    p = subprocess.Popen(f"module load {module} && echo $PATH", shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout, stderr = p.communicate()
    # the new PATH write into stdout
    path = stdout.decode()
    # print(old_env["PATH"])
    # print(current_env)
    env["PATH"] = path
    return env

def timeit(start):
    end = time.time()
    t = end-start
    if t < 60:
        return '{:.2f} seconds elapsed'.format(t)
    elif t < 3600:
        return '{:.2f} minutes elapsed'.format(t/60)
    else:
        return '{:.2f} hours elapsed'.format(t/3600)

def rank_worker(size, fill_data):
    worker_tasks = {w:[] for w in range(size)}
    w_idx = 0
    for fd in fill_data:
        worker_tasks[w_idx].append(fd)
        w_idx = (w_idx + 1) % size
    return worker_tasks

