# Sample specific weighted gene co-expression analysis for treatment outcome prediction in patients with metastatic clear cell renal cell carcinoma

## Abstract
Immunotherapies based on checkpoint blockade have recently become the standard of care for the treatment of metastatic clear cell renal cell carcinoma, the most prevalent subtype of kidney cancer. Although genomic alterations and immune microenvironments have been characterized in terms of their relationship with treatment response, gene network features remain absent in elucidating the driver and resistance to these therapies. We used transcriptomics data to conduct sample-specific weighted gene coexpression network inference and to explore the heterogeneity of patient clinical outcomes based on network features and biological pathway topology. Our results showed that network distance, adjusted with prior clinical knowledge, could measure the association of patients with treatment outcomes. Furthermore, higher gene connectivity and negative gene-gene associations were found to be related to poor prognosis and cancer recurrence post-therapy. In addition, entropy and topology scores calculated for biological pathways enabled the identification of the perturbed molecular processes associated with treatment responses that were not directly detected by gene expression changes. Overall, our sample-specific gene network analysis provides a comprehensive view on the capacity of network features for the stratification of cancer patients in the context of survival studies and prediction of response to treatment.


## SetPA (Sample specific entropy and topological score of Pathways)

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
