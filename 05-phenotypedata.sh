# Create covariate file (age, sex, predicted smoking, predicted cell counts
# Estimate first 10 principal components using GCTA

# split into males and females
# transform each sex using inverse normal
# recombine scaled to population mean and variance

# Decide on a format
# probably easiest to use plink format but check with Hash if this is good for matrixexqtl





# Normalise height and BMI

# Generate age acceleration residuals
Rscript resources/dnamage/dnamage.r ${beta_27k} ${phenfile} ${outfile}

# Generate smoking predictor
Rscript resources/smoking_predictor.R ${beta} ${phenfile} ${outfile}

# Generate cell counts
