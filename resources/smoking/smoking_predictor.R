# Thanks to Hannah Elliott for sending the original version of this script

arguments <- commandArgs(T)

methylationfile <- arguments[1]
phenfile <- arguments[2]
outfile <- arguments[3]

# Dataframe of Illig data (supplementary table 2) http://www.ncbi.nlm.nih.gov/pubmed/23691101
load("resources/illig.RData")

# Load methylation data - this contains the object "mbeta"
# This is an ncpg (rows) x nid (cols) matrix
# Columns and rows labelled with ID and CPG respectively
load(methylationfile)


predict.smoking <- function(Illig_data, mbeta)
{
	# Remove CpGs missing from data
	Illig_data <- subset(Illig_data, cpgs %in% rownames(mbeta))

	# Split data into CPG sites that increase or decrease in smokers
	Illig_data_up <- subset(Illig_data, Illig_data$all_effect >=0)
	Illig_data_down <- subset(Illig_data, Illig_data$all_effect <0)

	# subset SABRE data and order
	mbeta_down <- mbeta[rownames(mbeta) %in% Illig_data_down$cpgs,]
	mbeta_up <- mbeta[rownames(mbeta) %in% Illig_data_up$cpgs, ]

	# sort Illig data by Cpg name
	Illig_data_up <- Illig_data_up[order(Illig_data_up$cpgs),]
	Illig_data_down <- Illig_data_down[order(Illig_data_down$cpgs),]

	# sort SABRE data by Cpg name
	mbeta_up <- mbeta_up[order(rownames(mbeta_up)),]
	mbeta_down <- mbeta_down[order(rownames(mbeta_down)),]


	# Calculate diffs between SABRE beta values and the reference for each CpG site
	matrix_up <- mbeta_up - Illig_data_up$reference_never_median_beta_all
	matrix_down <- Illig_data_down$reference_never_median_beta_all - mbeta_down

	# Calculate scores
	scores_up <- c(t(matrix_up) %*% Illig_data_up$weights)
	scores_down <- c(t(matrix_down) %*% Illig_data_down$weights)

	# Combine scores
	scores_combined <- scores_up + scores_down

	dat <- data.frame(IID = rownames(mbeta), Smoking = scores_combined)
	return(dat)
}

smok <- predict.smoking(Illig_data, mbeta)

phen <- read.table(phenfile, header=T)
phen$index <- 1:nrow(phen)
phen <- merge(phen, smok, by="IID", all.x=T)
phen <- phen[order(phen$index),]
phen <- subset(phen, select=-c(index))

write.table(phen, file=outfile, row=F, col=T, qu=F)

