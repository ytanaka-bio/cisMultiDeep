import sys
import argparse
import warnings
warnings.filterwarnings("ignore")

#0.1 parse argument
parser = argparse.ArgumentParser(description='Identify the best feature set that can explain cell type',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
list_of_method = ['RFECVLinear','RFECVRandom','LassoCV','LassoLarsCV']
list_of_test = ["RandomForest","LinearRegress","skip"]

#0.2 subparameters mode
parser.add_argument('-g',nargs=1,help='how many top and bottom genes are used?',type=int,dest="num_gene",default=500)
parser.add_argument('-m',nargs=1,help='A method for feature selection',choices=list_of_method,dest="method",type=str,default=list_of_method[2])
parser.add_argument('-t',nargs=1,help='A method to assess feature selection, if you do not want to test, please choose "skip"',choices=list_of_test,dest="test",default=list_of_test[2])
parser.add_argument('-l',nargs=1,help='Lamda value for Lasso',type=float,dest="lamda",default=0.1)
parser.add_argument('-s',nargs=1,help='ratio of test dataset',default=0.2, type=float, dest="ratio_test")
parser.add_argument('-p',nargs=1,help='Number of threads',default=1,type=int,dest='thread')
parser.add_argument('-n',nargs="?",help='Do normalization?',default=True,dest='norm')
parser.add_argument('-f',nargs='+',help='H5AD file name',type=str,dest='h5ad_file',required=True)
parser.add_argument('-d',nargs=1,help='File name of differential feature analysis',type=str,dest="dif_file",required=True)
parser.add_argument('-a',nargs=1,help='Obs name including cell type information',type=str,dest="obs_name",required=True)
parser.add_argument('-r',nargs=1,help='Orthologous feature file name',type=str,dest='ortho',required=True)
parser.add_argument('-o',nargs=1,help='Prefix of output file name',type=str,dest="out_prefix",required=True)


#0.3 get arguments
args = parser.parse_args()
if args.method.__class__ != list:
    args.method = [args.method]
if args.test.__class__ != list:
    args.test = [args.test]

#0.4Import modules
import numpy as np
import pandas as pd
import statsmodels.api as sm
import scanpy as sc
from sklearn.feature_selection import RFECV
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LassoCV
from sklearn.linear_model import LassoLarsCV
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt


#definition
def model_assess(a, b, c, d, test):
    #Test
    if test == "RandomForest":
        model = RandomForestRegressor()       
    elif test == "LinearRegress":
        model = LinearRegression()
    else:
        sys.exit("Selected evaluation method is not supported\n");
    model.fit(a, b)
    r_sq = model.score(c,d)
    return r_sq
   
def feature_selection(a, b, c, d, method, test):
    if method == 'RFECVLinear':
        rfe = RFECV(LinearRegression(),n_jobs=args.thread[0])
        a_rfe = rfe.fit(a,b)
        selected = a.columns[a_rfe.support_]
        a_selected = a[selected]
        c_selected = c[selected]
    elif method == 'RFECVRandom':
        rfe = RFECV(RandomForestRegressor(),n_jobs=args.thread[0])
        a_rfe = rfe.fit(a,b)
        selected = a.columns[a_rfe.support_]
        a_selected = a[selected]
        c_selected = c[selected]
    elif method == 'LassoCV':
        ss = StandardScaler()
        a_std = ss.fit_transform(a)
        c_std = ss.transform(c)
        ls = LassoCV(n_jobs=args.thread[0],precompute=False)
        ls.fit(a_std,b)
        selected = a.columns[ls.coef_!=0]
        a_selected = a[selected]
        c_selected = c[selected]
    elif method == 'LassoLarsCV':
        ss = StandardScaler()
        a_std = ss.fit_transform(a)
        c_std = ss.transform(c)
        ls = LassoLarsCV(n_jobs=args.thread[0],precompute=False)
        ls.fit(a_std,b)
        selected = a.columns[ls.coef_!=0]
        a_selected = a[selected]
        c_selected = c[selected]
 
    else:
        sys.exit("Selected feature selection method is not supported\n");

    #print(a_selected)
    if test != "skip":
        r_sq_before = model_assess(a, b, c, d, test)
        r_sq_after  = model_assess(a_selected, b, c_selected, d, test)
    else:
        r_sq_before = 0
        r_sq_after = 0
        
    return a_selected, r_sq_before, r_sq_after
       
#read data files
df = pd.read_csv(args.dif_file[0],index_col=0)
celltype = df.columns
ortho = pd.read_csv(args.ortho[0])
ortho.columns = args.h5ad_file

combine = pd.DataFrame()
obs_data = pd.DataFrame()
for file in args.h5ad_file:
    print("Reading", file)
    adata = sc.read_h5ad(file)
    if args.norm == True:
        sc.pp.normalize_per_cell(adata,counts_per_cell_after=1e4)
        sc.pp.log1p(adata)

    #subset
    raw = adata[:,ortho[file]]
    data = raw.to_df()
    data = data.fillna(0)

    #combine
    data.columns = ortho[args.h5ad_file[0]]
    data = [combine,data]
    combine = pd.concat(data)

    data = [obs_data,raw.obs]
    obs_data = pd.concat(data)

#print(combine.shape)
#print(obs_data.shape)
Y = pd.DataFrame(0, index=combine.index,columns=celltype)

selected_features = []
for x in celltype:
    print("Analyzing",x,"by",*args.method)
    Y[x][obs_data[args.obs_name[0]] == x] = 1
    top_gene = df.sort_values(x,ascending=True).index[0:args.num_gene[0]]
    bottom_gene = df.sort_values(x,ascending=False).index[0:args.num_gene[0]]

    combine_top = combine[top_gene]
    combine_bottom = combine[bottom_gene]

    X = pd.concat([combine_top,combine_bottom],axis=1)
    X_train, X_test, Y_train, Y_test = train_test_split(X,Y[x],test_size=args.ratio_test,random_state=123)
    X_selected, r_sq_before, r_sq_after = feature_selection(X_train, Y_train, X_test, Y_test, *args.method, *args.test)
    if args.test[0] != "skip":
        print(x)
        print("R2 original       :", r_sq_before)
        print("R2 after selection:", r_sq_after)

    selected_features.extend(X_selected)

selected_features = list(set(selected_features))
combine_selected = combine[selected_features]

#output
output1 = args.out_prefix[0] + '_input.csv'
combine_selected.to_csv(output1)
output2 = args.out_prefix[0] + '_output.csv'
Y.to_csv(output2)
