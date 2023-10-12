# cisMultiDeep: Identifying Cell Type-Specific Cis-Regulatory Regions by Automatically-Tuned Deep Neural Network and SHAP

## Introduction
This repository presents the workflow to identify functional enhancers for each annotated cell type from a diverse collection of multi-omics profiles. Here, we hypothesize that cell type is determined by a common set of genes, whose expression is regulated by highly-conserved cis-regulatory regions across species. 

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
1.2. Add "chr" into Macaque bedpe files:
```{r eval=FALSE}
mkdir snm3C/Macaque2
mkdir snm3C/Macaque2/HiC_Loops
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Astro.loop.bedpe > snm3C/Macaque2/HiC_Loops/Astro.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/CLA.loop.bedpe > snm3C/Macaque2/HiC_Loops/CLA.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L23.loop.bedpe > snm3C/Macaque2/HiC_Loops/L23.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L4.loop.bedpe > snm3C/Macaque2/HiC_Loops/L4.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L5-ET.loop.bedpe > snm3C/Macaque2/HiC_Loops/L5-ET.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L5-IT.loop.bedpe > snm3C/Macaque2/HiC_Loops/L5-IT.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L6b.loop.bedpe > snm3C/Macaque2/HiC_Loops/L6b.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L6-CT.loop.bedpe > snm3C/Macaque2/HiC_Loops/L6-CT.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/L6-IT.loop.bedpe > snm3C/Macaque2/HiC_Loops/L6-IT.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Lamp5.loop.bedpe > snm3C/Macaque2/HiC_Loops/Lamp5.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/MG.loop.bedpe > snm3C/Macaque2/HiC_Loops/MG.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/NP.loop.bedpe > snm3C/Macaque2/HiC_Loops/NP.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/ODC.loop.bedpe > snm3C/Macaque2/HiC_Loops/ODC.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/OPC.loop.bedpe > snm3C/Macaque2/HiC_Loops/OPC.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Pvalb-BC.loop.bedpe > snm3C/Macaque2/HiC_Loops/Pvalb-BC.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Pvalb-ChC.loop.bedpe > snm3C/Macaque2/HiC_Loops/Pvalb-ChC.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Sncg.loop.bedpe > snm3C/Macaque2/HiC_Loops/Sncg.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Sst.loop.bedpe > snm3C/Macaque2/HiC_Loops/Sst.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Vip.loop.bedpe > snm3C/Macaque2/HiC_Loops/Vip.loop.bedpe 
awk 'OFS="\t" {$1="chr"$1; $4="chr"$4; print}' snm3C/Macaque/HiC_Loops/Vsc.loop.bedpe > snm3C/Macaque2/HiC_Loops/Vsc.loop.bedpe
```
1.3. Download gene coordinate GTF files as follow:
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
1.4. Remove duplicated genes from BED file:
```{r eval=FALSE}
R CMD BATCH remove_duplicate_gene_bed.R
```
1.5. Obtain orthologous gene list from [Biomart](http://useast.ensembl.org/biomart/martview) as follow:
- Choose "Ensembl Gene 10" and "Human genes (GRCh39.p14)"
- In Attributes section, choose "Homologues (Max select 6 orthologues)"
- In GENE tab, choose "Gene stable ID" and "Gene name"
- In ORTHOLOGUES [K-O] tab, choose "Macaque gene table ID", "Macaque gene name", "Mouse gene stable ID", and "Mouse gene name"
- In ORTHOLOGUES [U-Z] tab, choose "White-tufted-ear marmoset gene stable ID" and "White-tufted-ear marmoset gene name"
- Click "Result"
- Choose "compressed file (.gz)" in export all results, and click "GO"

### 2. Calculate the cell type specificity score for each gene and peak
2.1. Calculate the cell type specificity score for each gene from transcriptome (RNA):
```{r eval=FALSE}
python identify_celltype_gene.py
```
2.2. Calculate the cell type specificity score for each gene from methylome (mCG, mCH):
```{r eval=FALSE}
python identify_celltype_methyl.py
```
2.3. Calculate the cell type specificity score for each peak from chromatin accessibility profiles (ATAC):
```{r eval=FALSE}
python identify_celltype_chromatin.py
```
``
### 3. Identification of conserved peaks across species
3.1. Download LiftOver UCSC Chain files:
```{r eval=FALSE}
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToHg38.over.chain.gz       #Mouse vs Human
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToRheMac10.over.chain.gz   #Mouse vs Macaque
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/liftOver/mm10ToCalJac4.over.chain.gz    #Mouse vs Marmoset
```
3.2. Generate BED format files of ATAC peaks:
```{r eval=FALSE}
R CMD BATCH get_atac_bed.R
```
3.3. Run LiftOver command to identify orthologous regions of Mouse ATAC peaks in other species' genomes:
```{r eval=FALSE}
liftOver Mouse_atac.bed mm10ToHg38.over.chain.gz Mouse_atac_Human.bed Mouse_atac_Human_unmapped.bed
liftOver Mouse_atac.bed mm10ToRheMac10.over.chain.gz Mouse_atac_Macaque.bed Mouse_atac_Macaque_unmapped.bed
liftOver Mouse_atac.bed mm10ToCalJac4.over.chain.gz Mouse_atac_Marmoset.bed Mouse_atac_Marmoset_unmapped.bed
```
3.4. Run intersect command of BEDTools to identify ATAC peaks in each species corresponding to Mouse ATAC peaks:
```{r eval=FALSE}
bedtools intersect -a Mouse_atac_Human.bed -b Human_atac.bed -wa -wb -f 0.5  > Mouse_atac_Human_overlap.bed
bedtools intersect -a Mouse_atac_Macaque.bed -b Macaque_atac.bed -wa -wb -f 0.5  > Mouse_atac_Macaque_overlap.bed
bedtools intersect -a Mouse_atac_Marmoset.bed -b Marmoset_atac.bed -wa -wb -f 0.5  > Mouse_atac_Marmoset_overlap.bed
```
### 4. Get the list and dataset of orthologous genes and peaks for Deep Learning
4.1. Get conserved gene/peak list using a R script `get_cons_data.R`:
```{r eval=FALSE}
R CMD BATCH get_cons_data.R
```
4.2. Prepare input and output dataset for Deep Learning (Here, we focus on top 600 differential genes/peaks in each cell type):
```{r eval=FALSE}
python prepare_dataset.py -f 10XMultiome/Mouse/Mouse_rna.h5ad 10XMultiome/Human/Human_rna.h5ad 10XMultiome/Macaque/Macaque_rna.h5ad 10XMultiome/Marmoset/Marmoset_rna.h5ad -d all_rna_dif_cons.csv -a subclass_Bakken_2022 -r cons_rna_list.csv -o rna_600 -g 600
python prepare_dataset.py -f 10XMultiome/Mouse/Mouse_atac.h5ad 10XMultiome/Human/Human_atac.h5ad 10XMultiome/Macaque/Macaque_atac.h5ad 10XMultiome/Marmoset/Marmoset_atac.h5ad -d all_atac_dif_cons.csv -a subclass_Bakken_2022 -r cons_atac_list.csv -o atac_600 -g 600
python prepare_dataset.py -f snm3C/Mouse/Mouse_mCG_gene_fractions.h5ad snm3C/Human/Human_mCG_gene_fractions.h5ad snm3C/Macaque/Macaque_mCG_gene_fractions.h5ad snm3C/Marmoset/Marmoset_mCG_gene_fractions.h5ad -d all_mCG_dif_cons.csv -a subclass_Bakken_2022 -r cons_mCG_list.csv -o mCG_600 -g 600 -n False
python prepare_dataset.py -f snm3C/Mouse/Mouse_mCH_gene_fractions.h5ad snm3C/Human/Human_mCH_gene_fractions.h5ad snm3C/Macaque/Maacque_mCH_gene_fractions.h5ad snm3C/Marmoset/Marmoset_mCH_gene_fractions.h5ad -d all_mCH_dif_cons.csv -a subclass_Bakken_2022 -r cons_mCH_list.csv -o mCH_600 -g 600 -n False
```
### 5. Deep learning and SHAP value calculation
5.1. Train Deep learning model and calculate the contribution (SHAP value) of each gene/peak to the segregation of cell types:
```{r eval=FALSE}
python DeepSHAP.py -i rna_600_input.csv -p rna_600_output.csv -o rna_600_deep -t 12
python DeepSHAP.py -i atac_600_input.csv -p atac_600_output.csv -o atac_600_deep -t 12
python DeepSHAP.py -i mCG_600_input.csv -p mCG_600_output.csv -o mCG_600_deep -t 12
python DeepSHAP.py -i mCH_600_input.csv -p mCH_600_output.csv -o mCH_600_deep -t 12
```
### 6. Identify peaks and genes within the same HiC loop
6.1. Prepare gene and peak BED file for 600 differential genes/peaks:
```{r eval=FALSE}
R CMD BATCH prepare_bed.R
```
6.2. Identify genes and peaks within each HiC loop:
```{r eval=FALSE}
mkdir inLoop

#Human ATAC
bedtools pairtobed -a snm3C/Human/HiC_Loops/Astro.loop.bedpe -b Human_atac_600.bed > inLoop/Astro.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/CLA.loop.bedpe -b Human_atac_600.bed > inLoop/CLA.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L23.loop.bedpe -b Human_atac_600.bed > inLoop/L23.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-ET.loop.bedpe -b Human_atac_600.bed > inLoop/L5-ET.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-IT.loop.bedpe -b Human_atac_600.bed > inLoop/L5-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/NP.loop.bedpe -b Human_atac_600.bed > inLoop/NP.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-CT.loop.bedpe -b Human_atac_600.bed > inLoop/L6-CT.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-IT.loop.bedpe -b Human_atac_600.bed > inLoop/L6-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6b.loop.bedpe -b Human_atac_600.bed > inLoop/L6b.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Lamp5.loop.bedpe -b Human_atac_600.bed > inLoop/Lamp5.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/MG.loop.bedpe -b Human_atac_600.bed > inLoop/MG.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/ODC.loop.bedpe -b Human_atac_600.bed > inLoop/ODC.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/OPC.loop.bedpe -b Human_atac_600.bed > inLoop/OPC.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Pvalb-BC.loop.bedpe -b Human_atac_600.bed > inLoop/Pvalb-BC.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sncg.loop.bedpe -b Human_atac_600.bed > inLoop/Sncg.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sst.loop.bedpe -b Human_atac_600.bed > inLoop/Sst.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vip.loop.bedpe -b Human_atac_600.bed > inLoop/Vip.loop_atac.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vsc.loop.bedpe -b Human_atac_600.bed > inLoop/Vsc.loop_atac.txt

#Human RNA
bedtools pairtobed -a snm3C/Human/HiC_Loops/Astro.loop.bedpe -b Human_rna_600.bed > inLoop/Astro.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/CLA.loop.bedpe -b Human_rna_600.bed > inLoop/CLA.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L23.loop.bedpe -b Human_rna_600.bed > inLoop/L23.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-ET.loop.bedpe -b Human_rna_600.bed > inLoop/L5-ET.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-IT.loop.bedpe -b Human_rna_600.bed > inLoop/L5-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/NP.loop.bedpe -b Human_rna_600.bed > inLoop/NP.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-CT.loop.bedpe -b Human_rna_600.bed > inLoop/L6-CT.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-IT.loop.bedpe -b Human_rna_600.bed > inLoop/L6-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6b.loop.bedpe -b Human_rna_600.bed > inLoop/L6b.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Lamp5.loop.bedpe -b Human_rna_600.bed > inLoop/Lamp5.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/MG.loop.bedpe -b Human_rna_600.bed > inLoop/MG.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/ODC.loop.bedpe -b Human_rna_600.bed > inLoop/ODC.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/OPC.loop.bedpe -b Human_rna_600.bed > inLoop/OPC.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Pvalb-BC.loop.bedpe -b Human_rna_600.bed > inLoop/Pvalb-BC.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sncg.loop.bedpe -b Human_rna_600.bed > inLoop/Sncg.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sst.loop.bedpe -b Human_rna_600.bed > inLoop/Sst.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vip.loop.bedpe -b Human_rna_600.bed > inLoop/Vip.loop_rna.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vsc.loop.bedpe -b Human_rna_600.bed > inLoop/Vsc.loop_rna.txt

#Human mCG
bedtools pairtobed -a snm3C/Human/HiC_Loops/Astro.loop.bedpe -b Human_mCG_600.bed > inLoop/Astro.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/CLA.loop.bedpe -b Human_mCG_600.bed > inLoop/CLA.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L23.loop.bedpe -b Human_mCG_600.bed > inLoop/L23.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-ET.loop.bedpe -b Human_mCG_600.bed > inLoop/L5-ET.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-IT.loop.bedpe -b Human_mCG_600.bed > inLoop/L5-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/NP.loop.bedpe -b Human_mCG_600.bed > inLoop/NP.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-CT.loop.bedpe -b Human_mCG_600.bed > inLoop/L6-CT.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-IT.loop.bedpe -b Human_mCG_600.bed > inLoop/L6-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6b.loop.bedpe -b Human_mCG_600.bed > inLoop/L6b.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Lamp5.loop.bedpe -b Human_mCG_600.bed > inLoop/Lamp5.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/MG.loop.bedpe -b Human_mCG_600.bed > inLoop/MG.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/ODC.loop.bedpe -b Human_mCG_600.bed > inLoop/ODC.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/OPC.loop.bedpe -b Human_mCG_600.bed > inLoop/OPC.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Pvalb-BC.loop.bedpe -b Human_mCG_600.bed > inLoop/Pvalb-BC.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sncg.loop.bedpe -b Human_mCG_600.bed > inLoop/Sncg.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sst.loop.bedpe -b Human_mCG_600.bed > inLoop/Sst.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vip.loop.bedpe -b Human_mCG_600.bed > inLoop/Vip.loop_mCG.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vsc.loop.bedpe -b Human_mCG_600.bed > inLoop/Vsc.loop_mCG.txt

#Human mCH
bedtools pairtobed -a snm3C/Human/HiC_Loops/Astro.loop.bedpe -b Human_mCH_600.bed > inLoop/Astro.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/CLA.loop.bedpe -b Human_mCH_600.bed > inLoop/CLA.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L23.loop.bedpe -b Human_mCH_600.bed > inLoop/L23.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-ET.loop.bedpe -b Human_mCH_600.bed > inLoop/L5-ET.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L5-IT.loop.bedpe -b Human_mCH_600.bed > inLoop/L5-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/NP.loop.bedpe -b Human_mCH_600.bed > inLoop/NP.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-CT.loop.bedpe -b Human_mCH_600.bed > inLoop/L6-CT.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6-IT.loop.bedpe -b Human_mCH_600.bed > inLoop/L6-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/L6b.loop.bedpe -b Human_mCH_600.bed > inLoop/L6b.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Lamp5.loop.bedpe -b Human_mCH_600.bed > inLoop/Lamp5.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/MG.loop.bedpe -b Human_mCH_600.bed > inLoop/MG.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/ODC.loop.bedpe -b Human_mCH_600.bed > inLoop/ODC.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/OPC.loop.bedpe -b Human_mCH_600.bed > inLoop/OPC.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Pvalb-BC.loop.bedpe -b Human_mCH_600.bed > inLoop/Pvalb-BC.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sncg.loop.bedpe -b Human_mCH_600.bed > inLoop/Sncg.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Sst.loop.bedpe -b Human_mCH_600.bed > inLoop/Sst.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vip.loop.bedpe -b Human_mCH_600.bed > inLoop/Vip.loop_mCH.txt
bedtools pairtobed -a snm3C/Human/HiC_Loops/Vsc.loop.bedpe -b Human_mCH_600.bed > inLoop/Vsc.loop_mCH.txt

#Macaque ATAC
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Astro.loop.bedpe -b Macaque_atac_600.bed > inLoop/Astro.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/CLA.loop.bedpe -b Macaque_atac_600.bed > inLoop/CLA.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L23.loop.bedpe -b Macaque_atac_600.bed > inLoop/L23.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-ET.loop.bedpe -b Macaque_atac_600.bed > inLoop/L5-ET.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-IT.loop.bedpe -b Macaque_atac_600.bed > inLoop/L5-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/NP.loop.bedpe -b Macaque_atac_600.bed > inLoop/NP.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-CT.loop.bedpe -b Macaque_atac_600.bed > inLoop/L6-CT.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-IT.loop.bedpe -b Macaque_atac_600.bed > inLoop/L6-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6b.loop.bedpe -b Macaque_atac_600.bed > inLoop/L6b.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Lamp5.loop.bedpe -b Macaque_atac_600.bed > inLoop/Lamp5.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/MG.loop.bedpe -b Macaque_atac_600.bed > inLoop/MG.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/ODC.loop.bedpe -b Macaque_atac_600.bed > inLoop/ODC.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/OPC.loop.bedpe -b Macaque_atac_600.bed > inLoop/OPC.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Pvalb-BC.loop.bedpe -b Macaque_atac_600.bed > inLoop/Pvalb-BC.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sncg.loop.bedpe -b Macaque_atac_600.bed > inLoop/Sncg.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sst.loop.bedpe -b Macaque_atac_600.bed > inLoop/Sst.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vip.loop.bedpe -b Macaque_atac_600.bed > inLoop/Vip.loop_atac.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vsc.loop.bedpe -b Macaque_atac_600.bed > inLoop/Vsc.loop_atac.txt

#Macaque RNA
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Astro.loop.bedpe -b Macaque_rna_600.bed > inLoop/Astro.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/CLA.loop.bedpe -b Macaque_rna_600.bed > inLoop/CLA.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L23.loop.bedpe -b Macaque_rna_600.bed > inLoop/L23.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-ET.loop.bedpe -b Macaque_rna_600.bed > inLoop/L5-ET.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-IT.loop.bedpe -b Macaque_rna_600.bed > inLoop/L5-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/NP.loop.bedpe -b Macaque_rna_600.bed > inLoop/NP.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-CT.loop.bedpe -b Macaque_rna_600.bed > inLoop/L6-CT.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-IT.loop.bedpe -b Macaque_rna_600.bed > inLoop/L6-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6b.loop.bedpe -b Macaque_rna_600.bed > inLoop/L6b.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Lamp5.loop.bedpe -b Macaque_rna_600.bed > inLoop/Lamp5.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/MG.loop.bedpe -b Macaque_rna_600.bed > inLoop/MG.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/ODC.loop.bedpe -b Macaque_rna_600.bed > inLoop/ODC.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/OPC.loop.bedpe -b Macaque_rna_600.bed > inLoop/OPC.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Pvalb-BC.loop.bedpe -b Macaque_rna_600.bed > inLoop/Pvalb-BC.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sncg.loop.bedpe -b Macaque_rna_600.bed > inLoop/Sncg.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sst.loop.bedpe -b Macaque_rna_600.bed > inLoop/Sst.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vip.loop.bedpe -b Macaque_rna_600.bed > inLoop/Vip.loop_rna.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vsc.loop.bedpe -b Macaque_rna_600.bed > inLoop/Vsc.loop_rna.txt

#Macaque mCG
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Astro.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Astro.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/CLA.loop.bedpe -b Macaque_mCG_600.bed > inLoop/CLA.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L23.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L23.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-ET.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L5-ET.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-IT.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L5-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/NP.loop.bedpe -b Macaque_mCG_600.bed > inLoop/NP.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-CT.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L6-CT.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-IT.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L6-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6b.loop.bedpe -b Macaque_mCG_600.bed > inLoop/L6b.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Lamp5.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Lamp5.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/MG.loop.bedpe -b Macaque_mCG_600.bed > inLoop/MG.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/ODC.loop.bedpe -b Macaque_mCG_600.bed > inLoop/ODC.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/OPC.loop.bedpe -b Macaque_mCG_600.bed > inLoop/OPC.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Pvalb-BC.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Pvalb-BC.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sncg.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Sncg.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sst.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Sst.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vip.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Vip.loop_mCG.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vsc.loop.bedpe -b Macaque_mCG_600.bed > inLoop/Vsc.loop_mCG.txt

#Macaque mCH
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Astro.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Astro.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/CLA.loop.bedpe -b Macaque_mCH_600.bed > inLoop/CLA.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L23.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L23.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-ET.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L5-ET.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L5-IT.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L5-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/NP.loop.bedpe -b Macaque_mCH_600.bed > inLoop/NP.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-CT.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L6-CT.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6-IT.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L6-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/L6b.loop.bedpe -b Macaque_mCH_600.bed > inLoop/L6b.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Lamp5.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Lamp5.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/MG.loop.bedpe -b Macaque_mCH_600.bed > inLoop/MG.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/ODC.loop.bedpe -b Macaque_mCH_600.bed > inLoop/ODC.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/OPC.loop.bedpe -b Macaque_mCH_600.bed > inLoop/OPC.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Pvalb-BC.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Pvalb-BC.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sncg.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Sncg.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Sst.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Sst.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vip.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Vip.loop_mCH.txt
bedtools pairtobed -a snm3C/Macaque2/HiC_Loops/Vsc.loop.bedpe -b Macaque_mCH_600.bed > inLoop/Vsc.loop_mCH.txt

#Marmoset ATAC
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Astro.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Astro.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/CLA.loop.bedpe -b Marmoset_atac_600.bed > inLoop/CLA.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L23.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L23.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-ET.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L5-ET.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-IT.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L5-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/NP.loop.bedpe -b Marmoset_atac_600.bed > inLoop/NP.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-CT.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L6-CT.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-IT.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L6-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6b.loop.bedpe -b Marmoset_atac_600.bed > inLoop/L6b.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Lamp5.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Lamp5.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/MG.loop.bedpe -b Marmoset_atac_600.bed > inLoop/MG.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/ODC.loop.bedpe -b Marmoset_atac_600.bed > inLoop/ODC.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/OPC.loop.bedpe -b Marmoset_atac_600.bed > inLoop/OPC.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Pvalb-BC.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Pvalb-BC.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sncg.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Sncg.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sst.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Sst.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vip.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Vip.loop_atac.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vsc.loop.bedpe -b Marmoset_atac_600.bed > inLoop/Vsc.loop_atac.txt

#Marmoset RNA
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Astro.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Astro.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/CLA.loop.bedpe -b Marmoset_rna_600.bed > inLoop/CLA.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L23.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L23.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-ET.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L5-ET.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-IT.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L5-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/NP.loop.bedpe -b Marmoset_rna_600.bed > inLoop/NP.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-CT.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L6-CT.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-IT.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L6-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6b.loop.bedpe -b Marmoset_rna_600.bed > inLoop/L6b.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Lamp5.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Lamp5.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/MG.loop.bedpe -b Marmoset_rna_600.bed > inLoop/MG.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/ODC.loop.bedpe -b Marmoset_rna_600.bed > inLoop/ODC.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/OPC.loop.bedpe -b Marmoset_rna_600.bed > inLoop/OPC.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Pvalb-BC.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Pvalb-BC.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sncg.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Sncg.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sst.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Sst.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vip.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Vip.loop_rna.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vsc.loop.bedpe -b Marmoset_rna_600.bed > inLoop/Vsc.loop_rna.txt

#Marmoset mCG
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Astro.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Astro.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/CLA.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/CLA.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L23.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L23.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-ET.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L5-ET.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-IT.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L5-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/NP.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/NP.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-CT.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L6-CT.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-IT.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L6-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6b.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/L6b.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Lamp5.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Lamp5.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/MG.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/MG.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/ODC.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/ODC.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/OPC.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/OPC.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Pvalb-BC.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Pvalb-BC.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sncg.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Sncg.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sst.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Sst.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vip.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Vip.loop_mCG.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vsc.loop.bedpe -b Marmoset_mCG_600.bed > inLoop/Vsc.loop_mCG.txt

#Marmoset mCH
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Astro.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Astro.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/CLA.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/CLA.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L23.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L23.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-ET.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L5-ET.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L5-IT.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L5-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/NP.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/NP.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-CT.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L6-CT.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6-IT.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L6-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/L6b.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/L6b.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Lamp5.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Lamp5.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/MG.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/MG.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/ODC.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/ODC.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/OPC.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/OPC.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Pvalb-BC.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Pvalb-BC.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sncg.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Sncg.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Sst.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Sst.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vip.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Vip.loop_mCH.txt
bedtools pairtobed -a snm3C/Marmoset/HiC_Loops/Vsc.loop.bedpe -b Marmoset_mCH_600.bed > inLoop/Vsc.loop_mCH.txt

#Mouse ATAC
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Astro.loop.bedpe -b Mouse_atac_600.bed > inLoop/Astro.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/CLA.loop.bedpe -b Mouse_atac_600.bed > inLoop/CLA.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L23.loop.bedpe -b Mouse_atac_600.bed > inLoop/L23.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-ET.loop.bedpe -b Mouse_atac_600.bed > inLoop/L5-ET.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-IT.loop.bedpe -b Mouse_atac_600.bed > inLoop/L5-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/NP.loop.bedpe -b Mouse_atac_600.bed > inLoop/NP.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-CT.loop.bedpe -b Mouse_atac_600.bed > inLoop/L6-CT.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-IT.loop.bedpe -b Mouse_atac_600.bed > inLoop/L6-IT.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6b.loop.bedpe -b Mouse_atac_600.bed > inLoop/L6b.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Lamp5.loop.bedpe -b Mouse_atac_600.bed > inLoop/Lamp5.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/MG.loop.bedpe -b Mouse_atac_600.bed > inLoop/MG.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/ODC.loop.bedpe -b Mouse_atac_600.bed > inLoop/ODC.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/OPC.loop.bedpe -b Mouse_atac_600.bed > inLoop/OPC.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Pvalb-BC.loop.bedpe -b Mouse_atac_600.bed > inLoop/Pvalb-BC.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sncg.loop.bedpe -b Mouse_atac_600.bed > inLoop/Sncg.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sst.loop.bedpe -b Mouse_atac_600.bed > inLoop/Sst.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vip.loop.bedpe -b Mouse_atac_600.bed > inLoop/Vip.loop_atac.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vsc.loop.bedpe -b Mouse_atac_600.bed > inLoop/Vsc.loop_atac.txt

#Mouse RNA
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Astro.loop.bedpe -b Mouse_rna_600.bed > inLoop/Astro.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/CLA.loop.bedpe -b Mouse_rna_600.bed > inLoop/CLA.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L23.loop.bedpe -b Mouse_rna_600.bed > inLoop/L23.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-ET.loop.bedpe -b Mouse_rna_600.bed > inLoop/L5-ET.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-IT.loop.bedpe -b Mouse_rna_600.bed > inLoop/L5-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/NP.loop.bedpe -b Mouse_rna_600.bed > inLoop/NP.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-CT.loop.bedpe -b Mouse_rna_600.bed > inLoop/L6-CT.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-IT.loop.bedpe -b Mouse_rna_600.bed > inLoop/L6-IT.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6b.loop.bedpe -b Mouse_rna_600.bed > inLoop/L6b.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Lamp5.loop.bedpe -b Mouse_rna_600.bed > inLoop/Lamp5.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/MG.loop.bedpe -b Mouse_rna_600.bed > inLoop/MG.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/ODC.loop.bedpe -b Mouse_rna_600.bed > inLoop/ODC.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/OPC.loop.bedpe -b Mouse_rna_600.bed > inLoop/OPC.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Pvalb-BC.loop.bedpe -b Mouse_rna_600.bed > inLoop/Pvalb-BC.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sncg.loop.bedpe -b Mouse_rna_600.bed > inLoop/Sncg.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sst.loop.bedpe -b Mouse_rna_600.bed > inLoop/Sst.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vip.loop.bedpe -b Mouse_rna_600.bed > inLoop/Vip.loop_rna.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vsc.loop.bedpe -b Mouse_rna_600.bed > inLoop/Vsc.loop_rna.txt

#Mouse mCG
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Astro.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Astro.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/CLA.loop.bedpe -b Mouse_mCG_600.bed > inLoop/CLA.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L23.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L23.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-ET.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L5-ET.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-IT.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L5-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/NP.loop.bedpe -b Mouse_mCG_600.bed > inLoop/NP.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-CT.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L6-CT.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-IT.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L6-IT.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6b.loop.bedpe -b Mouse_mCG_600.bed > inLoop/L6b.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Lamp5.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Lamp5.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/MG.loop.bedpe -b Mouse_mCG_600.bed > inLoop/MG.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/ODC.loop.bedpe -b Mouse_mCG_600.bed > inLoop/ODC.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/OPC.loop.bedpe -b Mouse_mCG_600.bed > inLoop/OPC.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Pvalb-BC.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Pvalb-BC.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sncg.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Sncg.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sst.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Sst.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vip.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Vip.loop_mCG.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vsc.loop.bedpe -b Mouse_mCG_600.bed > inLoop/Vsc.loop_mCG.txt

#Mouse mCH
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Astro.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Astro.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/CLA.loop.bedpe -b Mouse_mCH_600.bed > inLoop/CLA.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L23.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L23.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-ET.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L5-ET.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L5-IT.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L5-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/NP.loop.bedpe -b Mouse_mCH_600.bed > inLoop/NP.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-CT.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L6-CT.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6-IT.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L6-IT.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/L6b.loop.bedpe -b Mouse_mCH_600.bed > inLoop/L6b.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Lamp5.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Lamp5.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/MG.loop.bedpe -b Mouse_mCH_600.bed > inLoop/MG.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/ODC.loop.bedpe -b Mouse_mCH_600.bed > inLoop/ODC.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/OPC.loop.bedpe -b Mouse_mCH_600.bed > inLoop/OPC.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Pvalb-BC.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Pvalb-BC.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sncg.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Sncg.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Sst.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Sst.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vip.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Vip.loop_mCH.txt
bedtools pairtobed -a snm3C/Mouse/HiC_Loops/Vsc.loop.bedpe -b Mouse_mCH_600.bed > inLoop/Vsc.loop_mCH.txt
```
## References
[BICCN Challenge](https://biccnchallenge.org/)
