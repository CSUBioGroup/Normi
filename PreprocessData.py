# -*- coding:utf-8 -*-
import pandas as pd

def write_file(filePath, text):
    file = open(filePath, 'w')
    sample = ''
    for t in text:
        sample = sample + str(t) + '\n'
    file.write(sample)
    file.close()


def smooth(df_p, df_e, slipe, k=5):
    '''
    Applying sliding size-fixed windows to expression profiles to obtain smoothed data

    Input: 
    df_p: dataframe of pseudo time 
    df_e: dataframe of expression profile with columns as genes, indexs as cells
    slipe: sliding length
    k: window size

    Output: a list containing the num of each linear trajectory and a dataframe of 
            expression profile after preprocessing
    '''
    count = []    
    for i in range(len(df_p.columns)):
        df_p_c = df_p.iloc[:, i:(i + 1)].dropna()
        df_p_c.columns = ['pseudotime']
        df_p_c.sort_values(by=['pseudotime'], ascending=True, inplace=True)
        df_p_e = pd.merge(df_p_c, df_e, left_index=True, right_index=True)
        df_ee = df_p_e.iloc[:, 1:]
        
        start = 0
        end = 0
        j = 0

        df_e_new = pd.DataFrame()
        while end==0:
            j = j + 1
            if (start+k) < len(df_ee.index):
                df = df_ee.iloc[start:(start+k), :]
            else:
                df = df_ee.iloc[(len(df_ee.index)-k):len(df_ee.index), :]
                end = 1

            # take zeros into consideration
            res = df.apply(lambda x:x.value_counts().get(0,0), axis=0).astype(float)
            for item in range(len(df.columns)):
                res[item]=0.0 if res[item]>=k/2 else df.iloc[:, item].mean()
            df_mean = res.to_frame().T
            df_mean.index=[j]

            if j==0:
                df_e_new = df_mean
            else:
                df_e_new = pd.concat([df_e_new, df_mean])

            start = start + slipe
        
        if i == 0:
            df_e_all = df_e_new
        else:
            df_e_all = pd.concat([df_e_all, df_e_new])

        count.append(len(df_e_new))
    return count, df_e_all