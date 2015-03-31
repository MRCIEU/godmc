# Create covariate file (age, sex, predicted smoking, predicted cell counts
# Estimate first 10 principal components using GCTA

# split into males and females
# transform each sex using inverse normal
# recombine scaled to population mean and variance

# Decide on a format
# probably easiest to use plink format but check with Hash if this is good for matrixexqtl


# Generate age acceleration residuals

Rscript resources/dnamage/dnamage.r ${beta_27k} resources/dnamage/probeAnnotation21kdatMethUsed.csv.gz resources/dnamage/datMiniAnnotation27k.csv.gz resources/dnamage/AdditionalFile3.csv.gz ${phenfile}




Save
 cell counts, batch effects (?), predicted smoking, predicted age, sex as covariate file in processed_data
 normalised height and bmi as phenotype file in processed_data 