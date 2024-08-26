import pandas as pd
import numpy as np 
from lifelines import CoxPHFitter
from contextlib import redirect_stdout
import io

"""
### generate common edges across samples
def get_sample_pair(sampleGroup):
    all_sample = []
    with open(f"../leave_one_out/{sampleGroup}.weight.txt") as f_w:
        for i in f_w:
            i= i.strip().split()
            all_sample.append(i[0])
    return all_sample

def get_common_edges(samGroup,sample_dir,cate):
    all_df = pd.read_csv(f"../leave_one_out/{sample_dir}/sample_specific_{sampleGroup[0]}.txt",\
                     sep="\t",index_col=None,header=None,names=["A","B",sampleGroup[0]])
    for i in sampleGroup[1:]:
        sample_df = pd.read_csv(f"../leave_one_out/{sample_dir}/sample_specific_{i}.txt",\
                            sep="\t",index_col=None,header=None,names=["A","B",i])
        all_df = pd.merge(all_df,sample_df,how="inner",on=["A","B"]) #outer"
        print(all_df.shape)
    all_df.to_csv(f"{cate}_above_common_edges.txt",index=None)

for i in ["primaryNivo","primaryEver","metaNivo","metaEver"]:
    sampleGroup = get_sample_pair(i)
    get_common_edges(sampleGroup,"all_sample_k_10_above",i)
"""

#start here to run analysis
def read_common_edges(cate):
    all_df= pd.read_csv(f"{cate}_common_edges.txt",index_col=None)
    df_cl = pd.read_csv("../leave_one_out/braun_data_clinical_data.csv")
    df_os = pd.merge(all_df.T.reset_index(),df_cl,how='left', left_on="index", right_on="RNA_ID")
    df_os_merge = df_os.T.copy()
    df_os_merge[1] = df_os_merge[0] +"_" + df_os_merge[1]
    df_os_merge = df_os_merge[df_os_merge.columns[1:]].T
    df_os_merge.columns = df_os_merge.iloc[0]
    df_os_merge.drop(1,inplace=True)
    df_os_merge.columns = df_os_merge.columns.tolist()[:-13] + df_os.columns.tolist()[-13:]
    df_os_merge.dropna(axis=1,inplace=True)
    
    for i in df_os_merge.iloc[:,1:-13].columns:
        up_fence  = df_os_merge[i].quantile(0.95)
        low_fence  = df_os_merge[i].quantile(0.05)
        df_os_merge[i][df_os_merge[i]>up_fence*1.1] = up_fence
        df_os_merge[i][df_os_merge[i]<low_fence*0.9] = low_fence
    return df_os_merge
    
def get_cox_result(df_os_input,cate):
    topvar = df_os_input.iloc[:,1:-13].var(axis=0).sort_values(ascending=False).index.tolist()[:50000]
    df_os_input.fillna(0,inplace=True)
    df_cox_os, df_cox_pfs = pd.DataFrame(), pd.DataFrame()
    df_os_input['OS_CNSR'].astype(bool)
    df_os_input['PFS_CNSR'].astype(bool)

    for i in topvar:
        try:
            df_os1 = df_os_input[[i,'OS', 'OS_CNSR']]
            cph = CoxPHFitter()
            cph.fit(df_os1, duration_col='OS', event_col='OS_CNSR')#,show_progress=True)
            df_cox_os = pd.concat([df_cox_os,cph.summary])

            df_os2 = df_os_input[[i,'PFS', 'PFS_CNSR']]
            cph = CoxPHFitter()
            cph.fit(df_os2, duration_col='PFS', event_col='PFS_CNSR')
            df_cox_pfs = pd.concat([df_cox_pfs,cph.summary]) 
        except:
            pass
    df_cox_os.to_csv(f"./cox_result/cox_result_{cate}_os1.txt")
    df_cox_pfs.to_csv(f"./cox_result/cox_result_{cate}_pfs1.txt")

for categroup in ["primaryNivo","primaryEver","metaNivo","metaEver"]: #
    df_os_temp = read_common_edges(categroup)
    #print(type(df_os_temp))
    get_cox_result(df_os_temp,categroup)
    