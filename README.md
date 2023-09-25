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
2. Obtain homologous gene list from [Biomart](http://useast.ensembl.org/biomart/martview) (save as `mart_export.txt.gz`)
- Choose "Ensembl Gene 10" and "Human genes (GRCh39.p14)
3. 
