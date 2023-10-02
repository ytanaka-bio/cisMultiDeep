# BICCN_Challenge_2023
## Introduction
This repository presents the workflow to identify functional enhancers for each annotated cell type from a diverse collection of multi-omics profiles. Here, we hypothesize that cell type is determined by a common set of genes, whose expression is regulated by highly-conserved enhancers across species. 

## Environment
```{r eval=FALSE}
python = 3.10.2
R = 4.0.6
```
## Requirement
### Python libraries
```{r eval=FALSE}
git clone https://github.com/ytanaka-bio/BICCN_Challenge_2023
cd BICCN_Challenge_2023
pip install -r requirements.txt --user
```
### R libraries
- [Seurat](https://satijalab.org/seurat/)
- [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html)
- [SpectralTAD](https://bioconductor.org/packages/release/bioc/html/SpectralTAD.html)
- [HiCcompare](https://www.bioconductor.org/packages/release/bioc/html/HiCcompare.html)
- [stringr](https://cran.r-project.org/web/packages/stringr/vignettes/stringr.html)

### Other softwares
- [AWS CLI](https://aws.amazon.com/jp/cli/)
- [Bedops](https://bedops.readthedocs.io/)
- [Bedtools](https://bedtools.readthedocs.io/en/latest/)
- [Kentutils](https://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads)

## Methods
### 1. Download and Preprocess datasets
1.1. Download cross-species' multi-omics datasets that were provided from BICCN committee using AWS CLI as follow:
```{r eval=FALSE}
aws s3 sync s3://biccn-challenge . --no-sign-request
```

1.2. Download gene coordinate GTF files as follow:
```{r eval=FALSE}
#Human
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz
zcat gencode.v33.annotation.gtf.gz | grep gene_name | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | sort -k1,1 -k2,2n -k3,3n > Human_gene.bed
mv gencode.v33.annotation.gtf.gz Human_gene.gtf.gz

#Macaque
wget https://ftp.ensembl.org/pub/release-110/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.110.gtf.gz
zcat Macaca_mulatta.Mmul_10.110.gtf.gz | sed -E 's/^(\w+)/chr\1/g' > Macaque_gene.gtf
rm Macaca_mulatta.Mmul_10.110.gtf.gz
grep gene_name Macaque_gene.gtf | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | sort -k1,1 -k2,2n -k3,3n > Macaque_gene.bed
gzip Macaque_gene.gtf 

#Marmoset
wget https://hgdownload.soe.ucsc.edu/goldenPath/calJac4/bigZips/genes/ncbiRefSeq.gtf.gz
zcat ncbiRefSeq.gtf.gz | grep gene_name | awk 'OFS="\t" {if ($3=="transcript") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | sort -k1,1 -k2,2n -k3,3n > Marmoset_gene.bed
mv ncbiRefSeq.gtf.gz Marmoset_gene.gtf.gz

#Mouse
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M22/gencode.vM22.annotation.gtf.gz
zcat gencode.vM22.annotation.gtf.gz | grep gene_name | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$14,$7}}' | tr -d '";' | sort -k1,1 -k2,2n -k3,3n > Mouse_gene.bed
mv gencode.vM22.annotation.gtf.gz Mouse_gene.gtf.gz
```

1.3. Obtain orthologous gene list from [Biomart](http://useast.ensembl.org/biomart/martview) as follow:
- Choose "Ensembl Gene 10" and "Human genes (GRCh39.p14)"
- In Attributes section, choose "Homologues (Max select 6 orthologues)"
- In GENE tab, choose "Gene stable ID" and "Gene name"
- In ORTHOLOGUES [K-O] tab, choose "Macaque gene table ID", "Macaque gene name", "Mouse gene stable ID", and "Mouse gene name"
- In ORTHOLOGUES [U-Z] tab, choose "White-tufted-ear marmoset gene stable ID" and "White-tufted-ear marmoset gene name"
- Click "Result"
- Choose "compressed file (.gz)" in export all results, and click "GO"

### 2. Calculate the cell type specificity score for each gene and peak
2.1. Calculate the cell type specificity score for each gene from transcriptome (RNA).
```{r eval=FALSE}
python identify_celltype_gene.py
```
2.2. Calculate the cell type specificity score for each gene from methylome (mCG, mCH).
```{r eval=FALSE}
python identify_celltype_methyl.py
```
2.3. Calculate the cell type specificity score for each peak from chromatin accessibility profiles (ATAC).
```{r eval=FALSE}
python identify_celltype_chromatin.py
```
### 3. Identification of conserved genes across species
```{r eval=FALSE}
R CMD BATCH get_cons_gene.R
```

### 4. Identification of converved peaks across species
4.1. Download LiftOver UCSC Chain files
```{r eval=FALSE}
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz       #Mouse vs Human
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToRheMac10.over.chain.gz   #Mouse vs Macaque
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToCalJac4.over.chain.gz    #Mouse vs Marmoset
```
4.2. Generate BED format files of ATAC peaks
```{r eval=FALSE}
R CMD BATCH get_atac_bed.R
```
4.3. Run LiftOver command to identify orthologous regions of Mouse ATAC peaks in other species' genomes.
```{r eval=FALSE}
liftOver Mouse_atac.bed mm10ToHg38.over.chain.gz Mouse_atac_Human.bed Mouse_atac_Human_unmapped.bed
liftOver Mouse_atac.bed mm10ToRheMac10.over.chain.gz Mouse_atac_Macaque.bed Mouse_atac_Macaque_unmapped.bed
liftOver Mouse_atac.bed mm10ToCalJac4.over.chain.gz Mouse_atac_Marmoset.bed Mouse_atac_Marmoset_unmapped.bed
```
4.4. Run intersect command of BEDTools to identify ATAC peaks in each species corresponding to Mouse ATAC peaks. 
```{r eval=FALSE}
bedtools intersect -a Mouse_atac_Human.bed -b Human_atac.bed -wa -wb -f 0.5  > Mouse_atac_Human_overlap.bed
bedtools intersect -a Mouse_atac_Macaque.bed -b Macaque_atac.bed -wa -wb -f 0.5  > Mouse_atac_Macaque_overlap.bed
bedtools intersect -a Mouse_atac_Marmoset.bed -b Marmoset_atac.bed -wa -wb -f 0.5  > Mouse_atac_Marmoset_overlap.bed
```
4.5. ee
## References
[BICCN Challenge](https://biccnchallenge.org/)
