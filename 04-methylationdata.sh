# Use R-GADA to generate structural variant data from methylation intensities
# Normalise methylation intensity data and adjust each probe for batch (fit as random effects or fixed effects ok?)
# For family data adjust methylation data for relatedness (take residuals after fitting pedigree matrix, i.e. GRAMMAR method)
# Use Houseman reference method to predict cell counts
# Use Zeilinger et al (2013) smoking associations to predict smoking status to be used as covariate in further analysis
# Create methylation relationship matrix for each chromosome and for all methylation probes in LDAK
# Adjust methylation data to remove outliers (rank normalisation)
# Filter methylation data on Naeem et al list (categories) and non-variant probes (threshold) and non-unique positions
# Create age predictor using 20 different methods provided by Steve Horvath


save
 structural variants in processed_data
 normalised methylation in processed_data
 
