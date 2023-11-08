import sys
import argparse
import warnings
warnings.filterwarnings("ignore")

#0.1 parse argument
parser = argparse.ArgumentParser(description='Prepare datasets for model training',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#0.2 subparameters mode
parser.add_argument('-g',help='how many top and bottom genes are used?',type=int,dest="num_gene",default=500)
parser.add_argument('-n',nargs="?",help='Do log normalization?',default=True,dest='norm')
parser.add_argument('-f',nargs=1,help='H5AD file name',type=str,dest='h5ad_file',required=True)
parser.add_argument('-d',nargs=1,help='File name of differential feature analysis',type=str,dest="dif_file",required=True)
parser.add_argument('-a',nargs=1,help='Obs name including cell type information',type=str,dest="obs_name",required=True)
parser.add_argument('-o',nargs=1,help='Prefix of output file name',type=str,dest="out_prefix",required=True)


#0.3 get arguments
args = parser.parse_args()

#0.4Import modules
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from scipy.sparse import csr_matrix
from scipy.sparse import save_npz

#read data files
df = pd.read_csv(args.dif_file[0],index_col=0)
celltype = df.columns

adata = sc.read_h5ad(args.h5ad_file[0])
if args.norm == True:
    sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
    sc.pp.log1p(adata)

#subset
raw = adata[:,df.index]
combine = raw.to_df()
combine = combine.fillna(0)

obs_data = raw.obs

#print(combine.shape)
#print(obs_data)
Y = pd.DataFrame(0, index=combine.index,columns=celltype)

selected_features = []
for x in celltype:
    print("Analyzing ", x)
    Y[x][obs_data[args.obs_name[0]] == x] = 1
    #bottom_gene = df.sort_values(x,ascending=True).index[0:args.num_gene]
    top_gene = df.sort_values(x,ascending=False).index[0:args.num_gene]
    selected_features.extend(top_gene)
    #selected_features.extend(bottom_gene)

selected_features = list(set(selected_features))
combine_selected = combine[selected_features]

#shape of input data
print(combine_selected.shape)

#output
output1 = args.out_prefix[0] + '_input.npz'
csr_matrix = csr_matrix(combine_selected.astype(pd.SparseDtype("float64",0)).sparse.to_coo())
save_npz(output1, csr_matrix)

output3 = args.out_prefix[0] + '_column.npy'
np.save(output3, combine_selected.columns)
output4 = args.out_prefix[0] + '_row.npy'
np.save(output4, combine_selected.index)

output2 = args.out_prefix[0] + '_output.csv'
Y.to_csv(output2)
