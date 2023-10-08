#module load python/3.10
import numpy as np
import pandas as pd
import scanpy as sc
import os
import matplotlib.pyplot as plt
from scipy.stats import ranksums

#make directory for differential methylation analysis
os.mkdir('methylation')

## mCG
Human_mCG_adata = sc.read_h5ad("snm3C/Human/Human_mCG_gene_fractions.h5ad")
raw_mCG = Human_mCG_adata.to_df()
celltype = Human_mCG_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Human_mCG_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Human_mCG_adata.obs.index[Human_mCG_adata.obs.subclass_Bakken_2022 == x]
    group2 = Human_mCG_adata.obs.index[Human_mCG_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCG.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCG.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Human_mCG_dif.csv")

Macaque_mCG_adata = sc.read_h5ad("snm3C/Macaque/Macaque_mCG_gene_fractions.h5ad")
raw_mCG = Macaque_mCG_adata.to_df()
celltype = Macaque_mCG_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Macaque_mCG_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Macaque_mCG_adata.obs.index[Macaque_mCG_adata.obs.subclass_Bakken_2022 == x]
    group2 = Macaque_mCG_adata.obs.index[Macaque_mCG_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCG.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCG.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Macaque_mCG_dif.csv")

Marmoset_mCG_adata = sc.read_h5ad("snm3C/Marmoset/Marmoset_mCG_gene_fractions.h5ad")
raw_mCG = Marmoset_mCG_adata.to_df()
celltype = Marmoset_mCG_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Marmoset_mCG_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Marmoset_mCG_adata.obs.index[Marmoset_mCG_adata.obs.subclass_Bakken_2022 == x]
    group2 = Marmoset_mCG_adata.obs.index[Marmoset_mCG_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCG.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCG.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Marmoset_mCG_dif.csv")

Mouse_mCG_adata = sc.read_h5ad("snm3C/Mouse/Mouse_mCG_gene_fractions.h5ad")
raw_mCG = Mouse_mCG_adata.to_df()
celltype = Mouse_mCG_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Mouse_mCG_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Mouse_mCG_adata.obs.index[Mouse_mCG_adata.obs.subclass_Bakken_2022 == x]
    group2 = Mouse_mCG_adata.obs.index[Mouse_mCG_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCG.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCG.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Mouse_mCG_dif.csv")

## mCH
Human_mCH_adata = sc.read_h5ad("snm3C/Human/Human_mCH_gene_fractions.h5ad")
raw_mCH = Human_mCH_adata.to_df()
celltype = Human_mCH_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Human_mCH_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Human_mCH_adata.obs.index[Human_mCH_adata.obs.subclass_Bakken_2022 == x]
    group2 = Human_mCH_adata.obs.index[Human_mCH_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCH.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCH.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Human_mCH_dif.csv")

Macaque_mCH_adata = sc.read_h5ad("snm3C/Macaque/Maacque_mCH_gene_fractions.h5ad")
raw_mCH = Macaque_mCH_adata.to_df()
celltype = Macaque_mCH_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Macaque_mCH_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Macaque_mCH_adata.obs.index[Macaque_mCH_adata.obs.subclass_Bakken_2022 == x]
    group2 = Macaque_mCH_adata.obs.index[Macaque_mCH_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCH.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCH.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Macaque_mCH_dif.csv")

Marmoset_mCH_adata = sc.read_h5ad("snm3C/Marmoset/Marmoset_mCH_gene_fractions.h5ad")
raw_mCH = Marmoset_mCH_adata.to_df()
celltype = Marmoset_mCH_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Marmoset_mCH_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Marmoset_mCH_adata.obs.index[Marmoset_mCH_adata.obs.subclass_Bakken_2022 == x]
    group2 = Marmoset_mCH_adata.obs.index[Marmoset_mCH_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCH.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCH.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Marmoset_mCH_dif.csv")

Mouse_mCH_adata = sc.read_h5ad("snm3C/Mouse/Mouse_mCH_gene_fractions.h5ad")
raw_mCH = Mouse_mCH_adata.to_df()
celltype = Mouse_mCH_adata.obs.subclass_Bakken_2022.unique()
score = pd.DataFrame(index=Mouse_mCH_adata.var_names,columns=celltype)

for x in celltype:
    group1 = Mouse_mCH_adata.obs.index[Mouse_mCH_adata.obs.subclass_Bakken_2022 == x]
    group2 = Mouse_mCH_adata.obs.index[Mouse_mCH_adata.obs.subclass_Bakken_2022 != x]
    def ranksums2(x):
        w = ranksums(x.loc[group1],x.loc[group2])
        return w
    def mean_dif(x):
        m = x.loc[group1].mean() - x.loc[group2].mean()
        return m
    test = raw_mCH.apply(ranksums2, axis=0)
    pval = -np.log10(test.loc[1,])
    ave  = raw_mCH.apply(mean_dif, axis=0)
    pval[ave < 0] = -pval[ave < 0]
    score[x] = pval

score.to_csv("methylation/Mouse_mCH_dif.csv")
