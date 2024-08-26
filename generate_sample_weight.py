#SIN weight by sample
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-i", type=str , default="./directory/expression.txt" , help="A path to 'gene expression matrix' file")
parser.add_argument("-o", type=str , default="./directory/sample_weight.txt" , help="A path to the output 'sample weight' file")

args = parser.parse_args()
input , output  = args.i , args.o

#datar = pd.DataFrame(np.triu(datar.values,1),index=datar.index,columns=datar.columns)
def get_weight(input, output):
    """
    generate sample weight (pearson correlation coefficient between the sample inside specified cohorts)
    """
    sample_name = input
    expr_data = pd.read_csv(sample_name, sep=',', index_col=0)
    datar = expr_data.T.corr(method='pearson',min_periods=10)
    data_stack = pd.DataFrame(datar.stack())
    data_stack.columns = ["ps"]
    data_stack = data_stack[data_stack.ps!= 0]
    data_stack = data_stack[data_stack.ps!= 1]
    r_max, r_min = np.max(data_stack.ps), np.min(data_stack.ps)
    sample_mean = data_stack.groupby(level=0).mean()
    sample_weight = (sample_mean -r_min  +0.01) / (r_max-r_min +0.01)
    sample_weight.to_csv(f"{output}",sep='\t',columns=None, index=True,header=False)

if __name__ == "__main__":
    get_weight(input, output)