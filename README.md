# BICCN_Challenge_2023
## Introduction
This repository presents the workflow to identify functional enhancers for each annotated cell type from a diverse collection of multi-omics profiles. Here, we hypothesize that cell type is determined by a common set of genes, whose expression is regulated by highly-conserved enhancers across species. 

## Environment
```{r eval=FALSE}
python = 3.10.2
R = 4.0.6
```
## Requirement
```{r eval=FALSE}
git clone https://github.com/ytanaka-bio/BICCN_Challenge_2023
cd BICCN_Challenge_2023
pip install -r requirements.txt --user
```
## Preprocessing datasets
1. Download cross-species' multi-omics datasets that were provided from BICCN committee using AWS CLI as follow:
```{r eval=FALSE}
aws s3 sync s3://biccn-challenge . --no-sign-request
```
2. Due to the limited number of HiC loops (<10,000 loops), four cell types (L5-ET, Pvalb-ChC, CLA, and Sncg) were removed from the subsequent analyses:
```{r eval=FALSE}
gzip snm3C/*/HiC_Loops/L5-ET.loop.bedpe
gzip snm3C/*/HiC_Loops/Pvalb-ChC.loop.bedpe
gzip snm3C/*/HiC_Loops/CLA.loop.bedpe
gzip snm3C/*/HiC_Loops/Sncg.loop.bedpe
```
3. Download genomic coordinate of genes with BED format using [UCSC Table Browser](https://genome.ucsc.edu/cgi-bin/hgTables) as follow:
- Human: GRCh38/hg38, All GENCODE V33, Basic, genome (region) (save as `hg38_gene.bed`)
- Macaque (Rhesus): Mmul_10/rheMac10, Ensembl Genes, ensGene, genome (region) (save as `rheMac10_gene.bed`) 
- Marmoset: Callithrix_jacchus_cj1700_1.1/calJac4, NCBI RefSeq, RefSeq All, genome (region) (save as `calJac4_gene.bed`)
- Mouse: GRCh38/mm10, All GENCODE VM22, Basic, genome (region) (save as `mm10_gene.bed`)

4. Obtain orthologous gene list from [Biomart](http://useast.ensembl.org/biomart/martview) as follow:
- Choose "Ensembl Gene 10" and "Human genes (GRCh39.p14)"
- In Attributes section, choose "Homologues (Max select 6 orthologues)"
- In GENE tab, choose "Gene stable ID" and "Gene name"
- In ORTHOLOGUES [K-O] tab, choose "Macaque gene table ID", "Macaque gene name", "Mouse gene stable ID", and "Mouse gene name"
- In ORTHOLOGUES [U-Z] tab, choose "White-tufted-ear marmoset gene stable ID" and "White-tufted-ear marmoset gene name"
- Click "Result"
- Choose "compressed file (.gz)" in export all results, and click "GO"
- Then, run a R script `get_cons_gene.R`,
```{r eval=FALSE}
R CMD BATCH get_cons_gene.R
```

## Calculate the cell type specificity score for each gene and peak
1. Calculate the cell type specificity score for each gene from transcriptome (RNA).
```{r eval=FALSE}
python identify_celltype_gene.py
```
2. Calculate the cell type specificity score for each gene from methylome (mCG, mCH).
```{r eval=FALSE}
python identify_celltype_methyl.py
```
3. Calculate the cell type specificity score for each peak from chromatin accessibility profiles (ATAC).
```{r eval=FALSE}
python identify_celltype_chromatin.py
```
4. ee
