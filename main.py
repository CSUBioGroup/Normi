 # -*- coding:utf-8 -*-
import pandas as pd
import argparse
from PreprocessData import *
from EstimateMI import *
from mRMR import *
from Evaluate import *


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--window_size', type=int, default=5, help='The size of sliding window to preprocess expressioon data.')
    parser.add_argument('--slide', type=int, default=1, help='Sliding length when preprocessing expression data.')
    parser.add_argument('--n_jobs', type=int, default=1, help='The num of available CPU to run this program.')
    parser.add_argument('--data_file', type=str, help='The input scRNA-seq gene expression file.')
    parser.add_argument('--time_file', type=str, help='The input time file.')
    parser.add_argument('--tf_file', type=str, default='', help='TF list.')
    parser.add_argument('--save_dir', type=str, default='./')

    opt = parser.parse_args()

    k = opt.window_size
    slide = opt.slide
    n_jobs = opt.n_jobs
    save_dir = opt.save_dir

    # --------------- STEP1 preprocess the dataset ---------------
    data_file = opt.data_file
    time_file = opt.time_file
    tf_file = opt.tf_file

    df_exp = pd.read_csv(data_file, header=0, index_col=0)
    df_pse = pd.read_csv(time_file, header=0, index_col=0)
    
    if tf_file != '':
        TF_set = set(pd.read_csv(tf_file, header=None, index_col=None).values.reshape(-1))
    else:
        TF_set = []

    count, df_exp = smooth(df_pse, df_exp, slide, k)

    # --------------- STEP2 estimate the MI ---------------
    df_mi = cal_mi2(df_exp, count, n_jobs, TF_set)
    df_mi = df_mi.loc[df_mi.score>0]
    df_mi.sort_values(by='score', ascending=False, inplace=True)

    # --------------- STEP3 filter out redundant edges using mRMR ---------------
    df_mrmr = MRMR2(df_mi, n_jobs)
    df_mrmr = df_mrmr.sort_values(by='score', ascending=False)
    df_mrmr = df_mrmr.loc[df_mrmr.score>0]
    df_mrmr.to_csv(save_dir + "/rankedEdges.csv", header=True, index=False)



if __name__ == '__main__':
    main()
