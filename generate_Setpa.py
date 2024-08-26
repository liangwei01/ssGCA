import argparse
import pandas as pd
import numpy as np 
import igraph as ig
import sys
from scipy.stats import zscore
import math

##############parser#############
parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-i", type=str , default="./directory/sample_weight.txt" , help="A path to 'sample weight' file")
parser.add_argument("-d", type=str , default="./ssGCN/" , help="A directory stores sample specific networks")
parser.add_argument("-g", type=str , default="c2.cp.kegg.v2023.1.Hs.symbols.gmt" , help="pathway gene sets from databsaes")
parser.add_argument("-score", type=str , default="entropy" , help="pathway scores: entropy or topology scores")


args = parser.parse_args()
sample_weight_file , sample_dir, gmt_file, pa_score = args.i , args.d, args.g, args.score

outfile = f"{pa_score}.txt"
#########################################



def read_g(sample_dir,sample):
    dftemp = pd.read_csv(f"{sample_dir}/sample_specific_{sample}.txt",sep="\t",header=None)
    dftemp.columns = [0,1,"weight"]
    dftemp["weight"] =abs(dftemp["weight"])
    all_genes = set(dftemp[0].tolist() + dftemp[1].tolist())
    sGCN = ig.Graph.TupleList(dftemp.itertuples(index=False),edge_attrs="weight")
    return sGCN, all_genes

def get_str(alist):
    return "_".join([str(i) for i in alist])


def get_sample_pair(sampleGroup):
    all_sample = []
    with open(f"{sampleGroup}") as f_w:
        for i in f_w:
            i= i.strip().split()
            all_sample.append(i[0])
    return all_sample

def trans_alist(alist):
    if type(alist) == int:
        return alist
    if type(alist) == float:
        return alist
    elif len(alist) > 0:
        thelist = pd.Series(alist)
        if thelist.sum() == 0:
            return 0
        else:
            #print(thelist.astype('float').sum(),thelist.count())
            return thelist.astype('float').sum()/thelist.count()
    else:
        return 0


def cal_entropy(alist):
    h = 0
    if sum(alist) != 0:
        temp = sum(alist)
    for i in alist:
        if i == 0:
            h += 0
        else:
            h += (i/temp) * math.log2(i/temp)
    return -h


if __name__ == "__main__":
    
    path_list = []
    with open(gmt_file) as gmt:
        for i in gmt:
            path_list.append([i.split()[0]] + i.split()[2:])

    sampleGroup = get_sample_pair(sample_weight_file)

    f_pathway = open(outfile ,"w")
    
    for i in sampleGroup:
        g = read_g(sample_dir,i)
        allg_genes = g[1]
        expected = []
        expected.append(i)
        for j in path_list:
            left_genes = set(j[1:]).intersection(allg_genes)
        ##extract subgraph based pathway genes
            subg = g[0].subgraph(left_genes)
            total, genes_nx = len(j[1:]), len(left_genes)
            #print([t for t in subg.es])
            if genes_nx > 1:
                if pa_score == "entropy":
                    eig = cal_entropy([es["weight"] for es in subg.es])
                elif pa_score == "eigenvector":
                    eig = subg.eigenvector_centrality(directed=False, scale=False, weights="weight")
                elif pa_score == "closeness":
                    eig =subg.closeness(mode='all', cutoff=None, weights="weight", normalized=True)
                elif pa_score == "edge_betweenness":
                    eig = subg.edge_betweenness(directed=False, cutoff=None, weights="weight")
                expected.append(str(trans_alist(eig)))
            else:
                expected.append(str(0))
        f_pathway.write("\t".join(expected)+"\n")
    f_pathway.close()