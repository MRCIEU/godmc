GoDMC pipeline
==============

# Data setup

## Genotype data

Imputed data to 1000 genomes version 1 phase 3
Alternative reference panels should be OK too.

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