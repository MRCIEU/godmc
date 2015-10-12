GoDMC pipeline
==============

# Data setup

## Genotype data

Imputed data to 1000 genomes version 1 phase 3
Alternative reference panels should be OK too. `rnorm(1000)`

Filtering:
- HWE < 1e-6
- MAF < 0.01
- Quality score < 0.8


cfgfg this is some code

    a <- rmorm(1000)

sdfsd


Format:
- Best guess genotypes (not dosages)
- Single binary plink format file (not separate chromosomes)


## Methylation data

Original IDAT files are required. They will then be normalised using R/meffil. 



## R packages

matrixeqtl
RPMM
plyr
dplyr
ggplot2
reshape2
knitr
rmarkdown
markdown
matrixStats
parallel
limma
nlme
quadprog
GenABEL

meffil maybe also





check data
- rename all SNPs to be chr_pos_a1_a2
- get allele frequencies of all SNPs

covariates
- create covariate dataset
- principal components
- cell counts
- smoking
- age
- sex


methylation
- normal betas
- rank transformed betas
- rank transformed and cell count adjusted betas
- methylation variance dataset

cnvs
- create cnv dataset










# Instructions

Clone the repository

	git clone git@scmv-ieugit.epi.bris.ac.uk:gh13047/godmc.git

Copy config.example to config and fill it in

run 01-runtest
 run 02-checkdata







add chr to snp list




1. rank transform predictor
2. lower and upper deciles
3. lower and upper quartiles
4. thresholds - below 5 vs above 5
5. thresholds - below 10 vs above 10
6. thresholds - below 5 vs about 10
7. k means into 2 groups


1. 1 for each cell type (9?) - rank transform after removing outsliers (± 5sd)
2. Diversity as predicted by shannon's entropy
