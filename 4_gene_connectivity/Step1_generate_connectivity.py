import pandas as pd
import sys

sampleGroup = sys.argv[1]
all_sample = []
with open(f"../leave_one_out/{sampleGroup}.weight.txt") as f_w:
    for i in f_w:
        i= i.strip().split()
        file = i[0]
        all_sample.append(file)

def get_connect(file_source):
    df_source = pd.read_csv(f"../leave_one_out/all_sample_k_10/sample_specific_{file_source}.txt",sep="\t",header=None)
    df_source.columns = ["node1", "node2","weight" ]
    df_source1 = df_source[["node2", "node1","weight" ]].copy()
    df_source1.columns = ["node1", "node2","weight" ]
    dfsum = pd.concat([df_source, df_source1]).drop_duplicates()
    df_source_table = dfsum.pivot(index='node1', columns='node2', values='weight')
    
    df_source_connectivity = df_source_table.fillna(0).sum().reset_index()
    df_source_connectivity.columns = ["Gene",file_source]
    return df_source_connectivity

df_all = get_connect(all_sample[0])

for i in all_sample[1:]:
    df_num = get_connect(i)
    df_all =pd.merge(df_all,df_num,on="Gene", how="outer")
    #print(df_all)

df_all.to_csv(f"{sampleGroup}.connectivity.txt",index=None)