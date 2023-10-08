#module load python/3.10
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt

#make directory for differential chromatin analysis
os.mkdir('chromatin_access')

#Differential expression across cell types in each species
Human_atac_adata = sc.read_h5ad("10XMultiome/Human/Human_atac.h5ad")
sc.pp.normalize_per_cell(Human_atac_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Human_atac_adata)
sc.tl.rank_genes_groups(Human_atac_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Human_atac_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Human_atac_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Human_atac_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]

score.to_csv('chromatin_access/Human_atac_dif.csv')
    
Macaque_atac_adata = sc.read_h5ad("10XMultiome/Macaque/Macaque_atac.h5ad")
sc.pp.normalize_per_cell(Macaque_atac_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Macaque_atac_adata)
sc.tl.rank_genes_groups(Macaque_atac_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Macaque_atac_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Macaque_atac_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Macaque_atac_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('chromatin_access/Macaque_atac_dif.csv')
 
Marmoset_atac_adata = sc.read_h5ad("10XMultiome/Marmoset/Marmoset_atac.h5ad")
sc.pp.normalize_per_cell(Marmoset_atac_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Marmoset_atac_adata)
sc.tl.rank_genes_groups(Marmoset_atac_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Marmoset_atac_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Marmoset_atac_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Marmoset_atac_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('chromatin_access/Marmoset_atac_dif.csv')

Mouse_atac_adata = sc.read_h5ad("10XMultiome/Mouse/Mouse_atac.h5ad")
sc.pp.normalize_per_cell(Mouse_atac_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Mouse_atac_adata)
sc.tl.rank_genes_groups(Mouse_atac_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Mouse_atac_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Mouse_atac_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Mouse_atac_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('chromatin_access/Mouse_atac_dif.csv')
 
