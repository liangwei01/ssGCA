# SetPA (Sample specific entropy and topological score of Pathways)

Sample specific entropy and topological scores for biological pathways

#Input File Format 
Gene expression matrix: columns-gene;  rows-patients

#Dependencies
The code is written in Python3. 
The packages are needed for the calculation: 
- Numpy 
- Pandas
- Igraph
- Scipy

#Basic Usage

Step1:
```
python generate_sample_weight.py -i motzer.log2.filter.txt -o motzer.log2.filter.weight.txt
```
`-i`: input file - gene expression matrix file
`-o`: output file - sample weight file

step2:
```
python construct.py -i motzer.log2.filter.txt -o ./ssGCN -s motzer.log2.filter.weight.txt
```

`-i`: input file - gene expression matrix file
`-s`: sample weight file
`-k`: balance paremeter from 0 to 1
`-o`: A path to store sample specific network matrix files

step3:
```
python generate_Setpa.py -i motzer.log2.filter.weight.txt -d ../sGCN354 -g c2.cp.kegg.v2023.1.Hs.symbols.gmt -score entropy
```
`-i`: input file - sample weight file
`-d`: a directory stores sample specific networks
`-g`: pathway gene sets from databsaes (https://www.gsea-msigdb.org/gsea/msigdb/)
`-score`: pathway scores: entropy or topology scores (entropy / eigenvector / closeness / edge_betweenness)
