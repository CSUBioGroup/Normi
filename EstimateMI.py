# -*- coding:utf-8 -*-
import numpy as np
import pandas as pd
import scipy.spatial as ss
from itertools import permutations, product
from pqdm.processes import pqdm
import dcor
from scipy.special import digamma
from math import log, floor


def MI_Gao(x,y,k=5):
    '''
    Estimate the mutual information I(X;Y) of X and Y from samples {x_i, y_i}_{i=1}^N
    Using *Mixed-KSG* mutual information estimator

    Input: 
    x: 2D array of size N*d_x (or 1D list of size N if d_x = 1)
    y: 2D array of size N*d_y (or 1D list of size N if d_y = 1)
    k: k-nearest neighbor parameter

    Output: one number of I(X;Y)
    '''

    assert len(x)==len(y), "Lists should have same length"
    assert k <= len(x)-1, "Set k smaller than num. samples - 1"
    N = len(x)
    data = np.c_[x, y]

    if x.ndim == 1:
        x = x.reshape((N,1))

    if y.ndim == 1:
        y = y.reshape((N,1))

    tree_xy = ss.cKDTree(data)
    tree_x = ss.cKDTree(x)
    tree_y = ss.cKDTree(y)

    knn_dis = [tree_xy.query(point,k+1,p=float('inf'))[0][k] for point in data]
    ans = 0

    for i in range(N):
        kp, nx, ny = k, k, k
        if knn_dis[i] == 0:
            kp = len(tree_xy.query_ball_point(data[i],1e-15,p=float('inf')))
            nx = len(tree_x.query_ball_point(x[i],1e-15,p=float('inf')))
            ny = len(tree_y.query_ball_point(y[i],1e-15,p=float('inf')))
            ans += log(kp) + log(N) - log(nx) - log(ny)
        else:
            nx = len(tree_x.query_ball_point(x[i],knn_dis[i]-1e-15,p=float('inf')))
            ny = len(tree_y.query_ball_point(y[i],knn_dis[i]-1e-15,p=float('inf')))
            ans += digamma(kp) + digamma(N) - digamma(nx) - digamma(ny)
    return ans/N


def cal_mi(i, j, x, y, count):
    '''
    Estimate the mutual information I(X;Y) of X and Y from samples {x_i, y_i}_{i=1}^N
    Using *Mixed-KSG* mutual information estimator

    Input: 
    i: name of TF
    j: name of TG
    x: the expression values of TF i
    y: the expression values of TG j
    count: a list including the num of cells in different branches
    MAXD: int, the maximum time lag

    Output: a dict with keys ['Gene1', 'Gene2']
    '''
    maxDelay = floor(min(count)/3)
    dc = []
    
    for d in range(maxDelay):
        xraw = x.copy()
        yraw = y.copy()
        sums = 0
        for k in range(len(count)):
            sums += count[k]

            # source gene
            idx_x = list(range(sums-(k+1)*maxDelay)) + list(range(sums-k*maxDelay, len(xraw)))
            xraw = xraw[idx_x]
            
            # target gene
            idx_y = list(range(sums-count[k]-k*maxDelay)) + list(range(sums-count[k]+d-k*maxDelay, len(yraw)))
            yraw = yraw[idx_y]
            idx_y = list(range(sums-maxDelay-k*maxDelay)) + list(range(sums-k*maxDelay-d, len(yraw)))
            yraw = yraw[idx_y]

        # calculate the dc, [0,1]
        d = dcor.distance_correlation(xraw, yraw)
        dc.append(d)
    dc = np.array(dc)
    MAXD = np.argmax(dc)  

    sums = 0
    # Xt_L:Xt-2, Xt-1, Xt
    # Yt_L:Yt-2, Yt-1, Yt
    for k in range(len(count)):
        xraw = x.copy()
        yraw = y.copy()
        for d in [0, MAXD]:
            idx = list(range(sums+d, sums+count[k]-MAXD+d))
            if d == 0:
                Xt_d = xraw[idx]
                if MAXD == 0:
                    Yt_d = yraw[idx]
                    break
            else:
                Yt_d = yraw[idx]
        if k == 0:
            Xt_L = Xt_d
            Yt_L = Yt_d
        else:
            Xt_L = np.r_[Xt_L, Xt_d]
            Yt_L = np.r_[Yt_L, Yt_d]
        sums += count[k]

    if Xt_L.ndim == 1:
        Xt_L = Xt_L.reshape((len(Xt_L),1))
    if Yt_L.ndim == 1:
        Yt_L = Yt_L.reshape((len(Yt_L),1))

    # I(Xt-L;Yt)
    mi = max(MI_Gao(Xt_L, Yt_L), 0)

    return {'Gene1':i, 'Gene2':j, 'score':mi}


def cal_mi2(data, count, n_jobs=1, TF_set=[]):
    
    print('---------- data.shape=', data.shape)

    if len(TF_set) == 0:
        gene_combs = list(permutations(data.columns.values, 2))
    else:
        TG_set = set(data.columns)
        gene_combs = product(TF_set, TG_set)
    gene_combs = filter(lambda x: x[0]!=x[1], gene_combs)
    params = [[v[0], v[1], data[v[0]].values, data[v[1]].values, count] for v in gene_combs]
    result = pqdm(params, cal_mi, n_jobs=n_jobs, argument_type='args', desc='Computations of MI')
    df_res = pd.DataFrame.from_dict(result)
    return df_res  
