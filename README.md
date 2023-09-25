# BICCN_Challenge_2023
## Introduction
Here, we present the workflow to identify functional enhancers for each annotated cell type from a diverse collection of multi-omics profiles.

## Environment
python = 3.10.2
R = 4.0.6

## Requirement
## Downloading
```{r eval=FALSE}
git clone https://github.com/ytanaka-bio/BICCN_Challenge_2023
cd BICCN_Challenge_2023
pip install -r requirements.txt --user
```
## Prepare datasets
1. Download multi-species' multi-omics datasets that were provided from BICCN committee using AWS CLI as follow:
 ```{r eval=FALSE}
aws s3 sync s3://biccn-challenge . --no-sign-request
```
2. Obtain orthologous gene list from [Biomart](http://useast.ensembl.org/biomart/martview) as follow: (save as `mart_export.txt.gz`)
- Choose "Ensembl Gene 10" and "Human genes (GRCh39.p14)
- In Attributes section, choose "Homologues (Max select 6 orthologues)"
- In GENE tab, choose "Gene stable ID" and "Gene name"
- In ORTHOLOGUES [K-O] tab, choose "Macaque gene table ID", "Macaque gene name", "Mouse gene stable ID", and "Mouse gene name"
- In ORTHOLOGUES [U-Z] tab, choose "White-tufted-ear marmoset gene stable ID" and "White-tufted-ear marmoset gene name"
- Click "Result"
- Choose "compressed file (.gz)" in export all results, and click "GO"
3. 
