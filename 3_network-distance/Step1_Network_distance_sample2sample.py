#!/usr/bin/env python
# coding: utf-8
import igraph as ig
import pandas as pd
from multiprocessing import Pool,get_context
from contextlib import redirect_stdout
import io

def read_g(file_name):
    df = pd.read_csv(file_name,sep="\t",header=None)
    df.columns = ["node1", "node2","weight" ]
    df["weight"] = abs(df["weight"])
    #df = df[df["weight"] > 0.7]
    GCN = ig.Graph.TupleList(df.itertuples(index=False),edge_attrs="weight")
    return GCN

def caculate_overlap(g1,g2):
    one_sample = read_g(f"./all_sample_k_10_above/sample_specific_{g1}.txt")
    another_sample = read_g(f"./all_sample_k_10_above/sample_specific_{g2}.txt")
    with redirect_stdout(io.StringIO()) as f:
        g_inter = ig.intersection([one_sample,another_sample]).ecount()
        g_union = one_sample.union(another_sample).ecount()
    return [g1, g2, str(g_inter/g_union)]

def get_sample_pair(sampleGroup):
    all_sample = []
    with open(f"{sampleGroup}.weight.txt") as f_w:
        for i in f_w:
            i= i.strip().split()
            all_sample.append(i[0])

    all_pair = []
    len_sample = len(all_sample)
    for i in range(len_sample):
        for j in range(i+1,len_sample):
            all_pair.append((all_sample[i],all_sample[j]))
    return all_pair

if __name__ == "__main__":
    sampleGroup = get_sample_pair("primaryNivo")
    f = open("network_distance_above_sample2sample.txt","w")
    with get_context("forkserver").Pool(25) as p:
        for result in p.starmap(caculate_overlap, sampleGroup):
            f.write(" ".join(result)+"\n")
