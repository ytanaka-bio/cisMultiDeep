#module load python/3.10
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt

#make directory for differential expression analysis
os.mkdir('gene_expression')

#Differential expression across cell types in each species
Human_rna_adata = sc.read_h5ad("10XMultiome/Human/Human_rna.h5ad")
sc.pp.normalize_per_cell(Human_rna_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Human_rna_adata)
sc.tl.rank_genes_groups(Human_rna_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Human_rna_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Human_rna_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Human_rna_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]

score.to_csv('gene_expression/Human_rna_dif.csv')
    
Macaque_rna_adata = sc.read_h5ad("10XMultiome/Macaque/Macaque_rna.h5ad")
sc.pp.normalize_per_cell(Macaque_rna_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Macaque_rna_adata)
sc.tl.rank_genes_groups(Macaque_rna_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Macaque_rna_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Macaque_rna_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Macaque_rna_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('gene_expression/Macaque_rna_dif.csv')
 
Marmoset_rna_adata = sc.read_h5ad("10XMultiome/Marmoset/Marmoset_rna.h5ad")
sc.pp.normalize_per_cell(Marmoset_rna_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Marmoset_rna_adata)
sc.tl.rank_genes_groups(Marmoset_rna_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Marmoset_rna_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Marmoset_rna_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Marmoset_rna_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('gene_expression/Marmoset_rna_dif.csv')

Mouse_rna_adata = sc.read_h5ad("10XMultiome/Mouse/Mouse_rna.h5ad")
sc.pp.normalize_per_cell(Mouse_rna_adata,counts_per_cell_after=1e4)
sc.pp.log1p(Mouse_rna_adata)
sc.tl.rank_genes_groups(Mouse_rna_adata, 'subclass_Bakken_2022', method='wilcoxon',key_added="wilcoxon",pts=True)
celltype = Mouse_rna_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Mouse_rna_adata.var_names,columns=celltype)
for x in celltype:
    test = sc.get.rank_genes_groups_df(Mouse_rna_adata,group=x,key="wilcoxon")
    test.index = test['names']
    pval = test['pvals']
    pval = -np.log10(pval)
    pval[test['logfoldchanges'] < 0] = -pval[test['logfoldchanges'] < 0]
    score[x] = pval[score.index]
    
score.to_csv('gene_expression/Mouse_rna_dif.csv')
 
