# -*- coding:utf-8 -*-
import pandas as pd
from sklearn import metrics


def concat_ref(df, df_r):
    '''
    Concat outputFile and referenceNetwork to get the ground truth of the edges

    Input: 
    df: dataFrame of the outputFile
    df_r: dataFrame of the reference network

    Output: dataFrame of the outputFile with ground truth
    '''
    df.columns = ['Gene1', 'Gene2', 'score']
    df_r.columns = ['Gene1', 'Gene2']
    df_r['real_score'] = 1
    df = pd.merge(df, df_r, how='left')
    df.loc[df.real_score.isna(), 'real_score'] = 0
    df = df.loc[df.Gene1!=df.Gene2].drop_duplicates(keep='first')
    return df


def cal_auc_aupr(df):
    '''
    Evaluate auroc and auprc of the algorithm on a specific dataset

    Input: 
    df: dataFrame of the outputFile with ground truth

    Output: a dict with keys ['auroc', 'auprc']
    '''
    scores = df.score.values
    true = df.real_score.astype(int).values
    fpr, tpr, thresholds = metrics.roc_curve(true, scores, pos_label=1)
    prec, recall, thresholds = metrics.precision_recall_curve(true, scores, pos_label=1)
    auroc = metrics.auc(fpr, tpr)
    auprc = metrics.auc(recall, prec)
    return {'AUROC':auroc, 'AUPRC':auprc}


def cal_EPR(df, k, sparsity):
    '''
    Evaluate EP, EPR of the algorithm on a specific dataset

    Input: 
    df: dataFrame of the inferred network with ground truth
    k: the top k edges are taken for evaluation 
    sparsity: the sparsity of the real network

    Output: a dict with keys of EP, EPR on the top k edges
    '''
    k = min(k, len(df))
    df = df.sort_values(by='score', ascending=False)
    TP = df.iloc[:k, :].real_score.sum().astype(int)
    EP = TP/k
    EPR = EP/sparsity
    return TP, EP, EPR