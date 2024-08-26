#!/usr/bin/env python
import sys
import pandas as pd
import numpy as np
import os

category, dir_k, k = sys.argv[1], sys.argv[2], float(sys.argv[3])

#expression_PxG_filename = "Expression_log2_21_" + category +"_braun_2020.csv"
expression_PxG_filename = "Expression_log2_21_" + category +"_braun_2020_above_1.csv"
# Load gene expression data
expr_data = pd.read_csv(expression_PxG_filename, sep=',', index_col=0)

#load cohort level network
datar = expr_data.corr(method='pearson',min_periods=1)


#Load Sample weights
def get_weight(NAME):
    dic_weight_TEMP = {}
    with open(f"{NAME}.above1.weight.txt") as f_w:
        for i in f_w:
            i= i.strip().split()
            dic_weight_TEMP[i[0]] = i[1]
    return dic_weight_TEMP

dic_weight= get_weight(category)


from scipy.stats import zscore
sample = list(expr_data.index )
sample_num = len(sample) -1

def filterByZscore(aNetwork_score):
    """
    apply z-score to filter weak edges
    """
    aNetwork_score = aNetwork_score.where(np.triu(np.ones(aNetwork_score.shape),1).astype(bool),np.nan)
    df_aNetwork = aNetwork_score.stack()
    df_aNetwork.index.names=['v1','v2']
    df_aNetwork = df_aNetwork.reset_index()
    df_aNetwork.columns=['feat1','feat2','weight']
    z_aNetwork = zscore(df_aNetwork["weight"])
    #GET R2 FILTERED BY Z-SCORE 2.58
    df_aNetwork.where(abs(z_aNetwork)>2.58,inplace=True)
    df_aNetwork = df_aNetwork.replace(0, np.nan)
    df_aNetwork.dropna(inplace=True)    
    return df_aNetwork

def generate_sample_edge(sample_p, all_sample_r):
    """
    generate sample specific network
    using pertuberd network - reference network
    the details can be accessed on sweet paper
    """
    #aggregate network
    expr_sample19 = expr_data.drop(sample_p,inplace=False)
    sample19 = expr_sample19.corr(method='pearson',min_periods=1)
    del expr_sample19
    #sample specific network
    sample_r = datar - sample19
    #sample specific confidence score
    #sample weight
    sam_weight = float(dic_weight[sample_p])
    sample_p_score = sam_weight * sample_num * k * (sample_r) + sample19
    del sample19 
    #sample_p = pd.DataFrame(np.triu(sample_p_score,1),index=datar.index,columns=datar.columns)
    #sample specific r2 > z-score
    df_sample_p = filterByZscore(sample_p_score)
    df_sample_p.to_csv(f"./{dir_k}/sample_specific_{sample_p}.txt",sep='\t',columns=None, index=False,header=False)

#generate cohort network
cohort_nn = filterByZscore(datar)
cohort_nn.to_csv(f"./{dir_k}/cohort_{category}.txt",sep='\t',columns=None, index=False,header=False)

#generate sample specific network
for i in sample:
    print(i)
    generate_sample_edge(i, datar)





