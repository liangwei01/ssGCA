#!/usr/bin/env python
# coding: utf-8
import igraph as ig
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import exp
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import time
from statistics import mean, median


def power_law(x1, a, b):
    return a*np.power(x1, b)

def get_r2(g):
    """
    generate the R2 rate for network scale free topology 
    """
    #number of edges: g.ecount()  number of vertices: g.vcount()
    print(g.ecount())
    sam_density = g.ecount()/(g.vcount()*(g.vcount()-1)/2)
    degree_sequence = sorted([d for d in g.degree()], reverse=True) 
    degree_out = {}
    for i in degree_sequence:
        if i in degree_out.keys():
            degree_out[i] += 1
        else:
            degree_out[i] = 1
    x = [int(i) for i in degree_out.keys()]
    y = [int(degree_out[i]) for i in x]
    popt, pcov = curve_fit(power_law, x, y)
    # get r2 via method 2
    r_squared_2 = r2_score(y, power_law(x, *popt), multioutput='variance_weighted')
    return [sam_density, r_squared_2]


##https://uppsala.instructure.com/courses/52162/pages/
##lab-network-construction-and-analysis-of-a-transcriptomic-and-metabolomic-dataset?module_item_id=310013
#function to get graph properties, takes a few minutes to run
def graph_prop(nn):
    """
    generate network features
    """
    a = time.time()
    ncount=nn.vcount()
    print(time.time()-a,"1 step")
    a = time.time()
    ecount=nn.ecount()
    print(time.time()-a,"2 step")
    a = time.time()
    kk = nn.degree()
    print(time.time()-a,"3 step")
    a = time.time()
    kk = [mean(kk),median(kk),max(kk)]
    print(time.time()-a,"4 step")
    a = time.time()
    """
    diameter=nn.diameter()
    print(time.time()-a,"5 step")
    a = time.time()
    av_path=nn.average_path_length()
    print(time.time()-a,"6 step")
    a = time.time()
    """
    dens=nn.density()
    print(time.time()-a,"7 step")
    """
    a = time.time()
    clustering=nn.transitivity_undirected() #this is the global clustering coefficient
    print(time.time()-a,"8 step")
    """
    a = time.time()
    conn=nn.is_connected()
    print(time.time()-a,"9 step")
    """
    a = time.time()
    min_cut=nn.mincut_value()
    print(time.time()-a,"10 step")
    """
    out = [ncount, ecount] + kk + [dens, conn]# clustering, min_cut
    #out=pd.DataFrame([ncount, ecount,kmean, kmedian,kmax, diameter, av_path, dens, clustering, conn, min_cut],
     #            index=['node_count','edge_count','kmean', 'kmedian','kmax','diameter','av_path_length','density','clustering_coef','connected?','minimum_cut']).T
    return out


def get_all_nn_feature(sample_file):
    df = pd.read_csv(sample_file,sep="\t",header=None)
    df.columns = ["node1", "node2","weight" ]
    df["weight"] = abs(df["weight"])
    sampleGCN = ig.Graph.TupleList(df.itertuples(index=False),edge_attrs="weight")
    all_result= get_r2(sampleGCN) + graph_prop(sampleGCN)
    return sample_file+"\t"+"\t".join([str(temp) for temp in all_result])+"\n"

def write_first_header(sampleGroup,dir_k):
    network_stats = open(f"{sampleGroup}_{dir_k}_network_stats.txt","w")
    headerline = ["sample_name","density", "r2",\
            'node_count','edge_count','kmean', \
              'kmedian','kmax',\
              'density','connected?']
    network_stats.write("\t".join(headerline)+"\n")
    network_stats.close()

# In[6]:


# multiproceesing must in    if __name__ == "__main__":
from multiprocessing import Pool, cpu_count
import sys 

if __name__ == "__main__":
    sampleGroup = sys.argv[1]
    dir_k = sys.argv[2]
    all_sample = []
    with open(f"{sampleGroup}.weight.txt") as f_w:
        for i in f_w:
            i= i.strip().split()
            file = f"./{dir_k}/sample_specific_{i[0]}.txt"
            all_sample.append(file)
    #add the cohort level GCN            
    all_sample.append(f"./{dir_k}/cohort_{sampleGroup}.txt")
    write_first_header(sampleGroup,dir_k)
    
    aaa = time.time()
    """
    #for loop 12 seconds for 3 samples
        result = []
    for  i in all_sample[:3]:
        result.append(get_all_nn_feature(i))
    """

    p = Pool(cpu_count()-1)
    result = []
    ###used method 4 seconds for 3 samples
    result = p.map(get_all_nn_feature, all_sample)
    p.close()
    p.join()
    network_stats = open(f"{sampleGroup}_{dir_k}_network_stats.txt","a")
    print([i for i in result])
    network_stats.write("".join(result))
    network_stats.close()
    print("used time: ", time.time()-aaa)
    print("finished")



