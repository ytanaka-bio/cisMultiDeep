# BICCN_Challenge_2023
## Introduction
This repository presents the workflow to identify functional enhancers for each annotated cell type from a diverse collection of multi-omics profiles. Here, we hypothesize that cell type is determined by a common set of genes, whose expression is regulated by highly-conserved enhancers across species. 

## Environment
```{r eval=FALSE}
python = 3.10.2
R = 4.0.6
```
## Requirement
###Python libraries
```{r eval=FALSE}
git clone https://github.com/ytanaka-bio/BICCN_Challenge_2023
cd BICCN_Challenge_2023
pip install -r requirements.txt --user
```
### R libraries
```{r eval=FALSE}
Seurat
Matrix
SpectralTAD
HiCcompare
stringr
```
### other softwares
```{r eval=FALSE}
bedops/2.4.41
```
## Downloading and Preprocessing datasets
1. Download cross-species' multi-omics datasets that were provided from BICCN committee using AWS CLI as follow:
```{r eval=FALSE}
aws s3 sync s3://biccn-challenge . --no-sign-request
```
2. Due to low quality of HiC data, four cell types (L5-ET, Pvalb-ChC, CLA, Sst, and Sncg) were removed from the subsequent analyses:
```{r eval=FALSE}
gzip snm3C/*/HiC_Loops/L5-ET.loop.bedpe
gzip snm3C/*/HiC_Loops/Pvalb-ChC.loop.bedpe
gzip snm3C/*/HiC_Loops/CLA.loop.bedpe
gzip snm3C/*/HiC_Loops/Sncg.loop.bedpe
gzip snm3C/*/HiC_Loops/Sst.loop.bedpe
```
3. Download genomic coordinate of genes with GTF format as follow:
```{r eval=FALSE}
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz
wget https://ftp.ensembl.org/pub/release-110/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.110.gtf.gz
zcat Macaca_mulatta.Mmul_10.110.gtf.gz | sed -E 's/^(\w+)/chr\1"/g' > Macaque_gene.gtf
rm Macaca_mulatta.Mmul_10.110.gtf.gz
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009663435.1/download?include_annotation_type=GENOME_GTF&filename=GCF_009663435.1.zip" -H "Accept: application/zip"
unzip GCF_009663435.1.zip 
sed -E 's/gene /gene_name /g' < ncbi_dataset/data/GCF_009663435.1/genomic.gtf > Marmoset_gene.gtf
rm GCF_009663435.1.zip
rm -rf ncbi_dataset/
```

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
