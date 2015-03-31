################################################
# SAMPLE CODE FOR HOUSEMAN, ACCOMANDO, ET AL.
# BMC Bioinformatics (2012)
#
#
# Version 1.10
# June 8, 2012
#
# Author: E. Andres Houseman
#
################################################

# Define wbcInference functions and load required libraries
source("resources/houseman/wbcInference-V112.R")
library(nlme)

# Load sample data [*Footnote 1*]
load("resources/houseman/Houseman-WBC-Data-V110.Rdata")

load(methylation_data)

phen <- data.frame(1, ncol(mbeta)))
rownames(phen) <- colnames(mbeta)


##############################################
##### ##### Step 1: Fit Validation Model (S0)
##### #####        (Required for all examples)
##############################################

validationData_Assay <- validationData_Assay[rownames(validationData_Assay) %in% rownames(mbeta), ]

validEst = validationWBC(
              validationData_Assay,             # Validation methylation (CpGs x subjects)
              validationData_Pheno,             # Validation phenotype frame (subjects x covariates)
              y~PanTcell+CD8T+CD4T+NK+Bcell+Mono+Gran, # Validation model (fixed effects) [*Footnote 2*]
              ~1|BeadChip                       # Validation batch adjustment (random effects) [*Footnote 3*]
           )
validEst

##################################################################################
##### ##### Example 1: Fit regression model to Target (S1)
##################################################################################

# Choose top pseudo-DMRs
DMRselection = rownames(validEst$coefEsts[validEst$orderFstat[1:100],])

# Check to make sure CpG names match!
all(rownames(validEst$coefEsts[DMRselection,]) == rownames(mbeta[DMRselection,]))

# Reporter matrix for PanT rearrangement [*Footnote 4*]
Lwbc = diag(8)[-(1:2),] 
Lwbc[1:2,2] = 1 ; Lwbc[,1] = 1
rownames(Lwbc) = colnames(validEst$coefEsts)[-(1:2)]
colnames(Lwbc) = colnames(validEst$coefEsts)

Lwbc # View reporter matrix

targetEst = inferWBCbyLme(
  mbeta[DMRselection,],     # Target methylation (CpGs x subjects)
  phen,               # Target phenotype frame (subjects x covariates)
  y~case+gender+ageCtr,                     # Target model (fixed effects) [*Footnote 2*]
  ~1|BeadChip,                              # Target batch adjustment (random effects) [*Footnote 3*]
  validEst$coefEsts[DMRselection,],         # Raw coefficient estimates for WBC 
  Lwbc                                      # Reporter matrix [*Footnote 4*]
)

targetEst # View model estimates

### Get bootstraps [*Footnote 3*]
# Warning:  this can take a long time
targetBoot = bootInferWBCbyLme( 
  mbeta[DMRselection,],     # Target methylation (CpGs x subjects)
  targetDataHNSCC_Covariates,               # Target phenotype frame (subjects x covariates)
  y~case+gender+ageCtr,                     # Target model (fixed effects) [*Footnote 2*]
  ~1|BeadChip,                              # Target adjustment (random effects) [*Footnote 5*] 
  validEst$coefEsts[DMRselection,],         # Raw coefficient estimates for WBC 
  Lwbc,                                     # Reporter matrix [*Footnote 4*]
  R=250,                                    # Number of bootstrap samples to run
  vcovWBC=validEst$coefVcovs[validEst$orderFstat[1:100]],   # WBC fixed effects v-cov estimates [*Footnote 6*]
  degFree=validEst$degFree[validEst$orderFstat[1:100]]      # WBC degrees-of-freedom [*Footnote 7*]
)

targetBoot # View bootstrap summary

# View summary with bootstraps
summary(targetEst, targetBoot)


##################################################################################
##### ##### Example 2: Conduct a permutation test
##################################################################################

########################################################
# Get regression estimate of inferred WBC on covariates
########################################################

OmegaEst = projectWBC( #Impute WBC by projection
  mbeta[DMRselection,],
  validEst$coefEsts[DMRselection,],    
  Lwbc)

GammaLM = lm(OmegaEst~case+gender+ageCtr, targetDataHNSCC_Covariates)
GammaCoef = GammaLM$coef
GammaPvals = sapply(summary(GammaLM), function(u)u$coef[,4])

GammaPvalCase0 = GammaPvals["case",]

# View estimates
# Note that these will not exactly match Example 1 because
# Example 1 takes into account batch effects.
round(100*GammaCoef,1)

########################################################
# Use permutation test (bootstrapping errors in validation)
#   to test association of case with WBC
########################################################

# Initialize parametric bootstrap for validation data
initParamBoot = initializeBootSampleValidation(validEst, validEst$orderFstat[1:100])

NPERMS = 1000
GammaPvalCasePerms = matrix(NA, NPERMS, length(GammaPvalCase0))
for(r in 1:NPERMS){
  OmegaBoot = projectWBC( #Impute WBC by projection
    mbeta[DMRselection,],
    bootSampleValidation(initParamBoot), # Parm bootstrap for validation data
    Lwbc)
  phenoPerm = targetDataHNSCC_Covariates
  phenoPerm$case = sample(phenoPerm$case)
  GammaPvalCasePerms[r,] = sapply(summary(lm(OmegaBoot~case+gender+ageCtr, phenoPerm)), function(u)u$coef["case",4]) 
}

# Omnibus (permutation) p-value
mean(min(GammaPvalCase0)>=apply(GammaPvalCasePerms,1,min))

##################################################################################
##### ##### Example 3: Project Target (S1) onto Validation (S0) and infer WBC
##################################################################################

####### Projections for HNSCC data

unconstrainedCoefs = projectWBC(
  mbeta[DMRselection,],
  validEst$coefEsts[DMRselection,],    
  Lwbc, nonnegative = FALSE)

constrainedCoefs = projectWBC(
  mbeta[DMRselection,],
  validEst$coefEsts[DMRselection,],    
  Lwbc)

head(unconstrainedCoefs)
head(constrainedCoefs)

# View heatmap of coefficients
heatmap(log(1+constrainedCoefs), scale="n", 
  col=hsv(0.75,seq(0,1,0.01),1),
  RowSide=c("green","red")[targetDataHNSCC_Covariates$case+1])

####### Projections for mixture experiment data

ObsMix = projectWBC(
  mixtureExperiment_Assay[DMRselection,],
  validEst$coefEsts[DMRselection,],    
  Lwbc)

ExMix = matrix(0, 12, 5)
colnames(ExMix) = c("Pan-T", colnames(ObsMix)[-(1:2)])
for(i in 1:12){
  ExMix[i,"Bcell"] = mixtureExperiment_Design$B[i]
  ExMix[i,"Pan-T"] = mixtureExperiment_Design$T[i]
  ExMix[i,"Gran"] = mixtureExperiment_Design$Gran[i]
  ExMix[i,"Mono"] = mixtureExperiment_Design$Mono[i]
}

colnames(ObsMix) = c("T (CD8+)", "T (CD4+)", "NK", "B Cell", "Monocyte", "Granulocyte")
colnames(ExMix) = c("T Cell", "NK", "B Cell", "Monocyte", "Granulocyte")
rownames(ObsMix) = rownames(ExMix) = mixtureExperiment_Design$strip

colSortO = c(5,4,3,1,2,6)[6:1]
colSortE = c(4,3,2,1,5)[5:1]
rowSort = c(2*(1:6)-1, 2*(1:6))

# Print Observed and Expected
round(100*ObsMix[rowSort, colSortO],1)
round(100*ExMix[rowSort, colSortE],1)

# Print errors
ObsMixCollapse = ObsMix[rowSort, colSortO]
ObsMixCollapse = cbind(ObsMixCollapse[,1], 
   ObsMixCollapse[,2]+ObsMixCollapse[,3], ObsMixCollapse[,4:6])
colnames(ObsMixCollapse) = colnames(ExMix)
round(100*(ObsMixCollapse-ExMix[rowSort, colSortE]),1)


############################################

##### ##### FOOTNOTES
#
#[*Footnote 1*] 
#
# Test data consist of 500 CpGs chosen from Illumina Infinium 27K array
# These were the 500 most informative CpGs for distinguishing WBC type
# (obtained by ordering the ANOVA F-statistics from largest to smallest
#  in the manner demonstrated in this file). Note that a LME is used to
#  adjust for batch effect, not ComBat, because LME error properties are
#  better understood.
#
# The CpGs are ordered from "most informative" to "least informative".
#
#
#[*Footnote 2*]
#
# Proper syntax for specifying a fixed effects model is as usual, except that
# the outcome must be specified as "y".  Consequently, "y" is a reserved name
# in this function (i.e. any variable in the phenotype file named "y" will
# not be accessible). 
#
#
#[*Footnote 3*]
#
# Batch adjustment by LME is obtained by specifying an appropriate
#  non-hierarchical random effects model definition (passed to nlme:::lme).
#  This will be slow for dense arrays. OLS will be used if this is omitted.
#  ComBat may be supported in the future, but not now. The proper syntax is 
#  "~1|BeadChip", where "BeadChip" is the name of the variable in the 
#  phenotype file that specifies chip number.
#
#
#[*Footnote 4*]
#
# The hierarchical design reflected in the WBC phenotype file
# requires that an appropriate "reporter" matrix be defined to extract desired 
# WBC types.  
#
#
#[*Footnote 5*]
#
# In this WBC function, LME-adjustment for batch effects is *required*
# Another function implements a faster, OLS-based algorithm that ignores
# batch effects ("inferWBCbyLm").  The proper syntax is "~1|BeadChip",
# where "BeadChip" is the name of the variable in the phenotype file that
# specifies chip number.
#
#
#[*Footnote 6*]
#
# The argument "vcovWBC" may be omitted if no double-bootstrap estimates
# are desired.  However, the double-bootstrap estimates are not particularly
# processor-intensive, so there is little additional cost to computing them.
# They do require WBC fixed effects variance-covariance estimates.
# Make sure to select the sub-list that corresponds to the CpGs selected for
# the computation!
#
#
#[*Footnote 7*]
#
# The argument "degFree" may be omitted if Gaussian noise is to be used
# in the double-bootstrap, instead of t-distributed noise.  T-distributed noise
# is likely to be (slightly) more accurate.
#
