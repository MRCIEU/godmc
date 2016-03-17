# Thanks to Hannah Elliott for sending the original version of this script


main <- function()
{
	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	fam_file <- arguments[2]
	out_file <- arguments[3]
	smoking.prediction.plot <- arguments[4]
	SD <- as.numeric(arguments[5])
	cov_file <- arguments[6]
	gwas_list_file <- arguments[7]
	
	
	# Dataframe of Illig data (supplementary table 2) http://www.ncbi.nlm.nih.gov/pubmed/23691101
	load("resources/smoking/illig.RData")
	fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]
	
	covs <- read.table(cov_file, header=T, stringsAsFactors=FALSE)
	m <- match(fam[,2],covs$IID)
	covs <- covs[m,]
	
	# Load methylation data - this contains the object "norm.beta"
	# This is an ncpg (rows) x nid (cols) matrix
	# Columns and rows labelled with ID and CPG respectively

	message("Loading methylation data")
	load(methylationfile)
	m <- match(fam[,2],colnames(norm.beta))
	norm.beta <- norm.beta[,m]

	message("Predicting smoking levels")
	smok <- predict.smoking(Illig_data, norm.beta)
	nom <- names(smok)
	smok <- merge(smok, fam, by.x="IID", by.y="V2")
	smok <- subset(smok, select=c("V1", nom))
	m <- match(fam$V2,smok$IID)
	smok <- smok[m,]

	pdf(smoking.prediction.plot, width=12, height=8)
	par(mfrow=c(2,2))
	plot(smok$Smoking, xlab="", main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),")",sep=""),cex.main=0.7)
	hist(smok$Smoking, xlab="", main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),")",sep=""),cex.main=0.7)
	abline(v=mean(smok$Smoking,na.rm=T)-SD*sd(smok$Smoking,na.rm=T),lty=2)
	abline(v=mean(smok$Smoking,na.rm=T)+SD*sd(smok$Smoking,na.rm=T),lty=2)
	qqnorm(smok$Smoking, main=paste("Smoking prediction (N=", length(which(!is.na(smok$Smoking))),"; shapiroP=",signif(as.numeric(shapiro.test(smok$Smoking)[2]),2),")",sep=""),cex.main=0.7)
	qqline(smok$Smoking)
	par(mfrow=c(2,2))
	plot(covs$Age_numeric,smok$Smoking, xlab="Age", ylab="predicted smoking",main="Age vs Smoking prediction",cex.main=0.7)
	
	quiet <- dev.off()

	# Rank transform
	smok$Smoking <- rntransform(smok$Smoking)

	write.table(subset(smok, select=-c(V1)), file=paste0(out_file, ".txt"), row=F, col=T, qu=F)
	write.table(smok, file=paste0(out_file, ".plink"), row=F, col=F, qu=F)

	# Get individuals with age >= 25

	gwas_list <- basename(out_file)

	i <- subset(covs, Age_numeric >= 25)$IID
	smoki <- subset(smok, IID %in% i)

	if(nrow(smoki) > 50)
	{
		gwas_list <- c(gwas_list, paste0(basename(out_file), "_ge25"))
		m <- "GWAS will be performed"
	} else {
		m <- "GWAS won't be performed"
	}

	message(nrow(smoki), " individuals over age 25 with predicted smoking data, ", m)
	write.table(smoki, file=paste0(out_file, "_ge25.plink"), row=F, col=F, qu=F)

	j <- subset(covs, Age_numeric < 25)$IID
	smokj <- subset(smok, IID %in% j)

	if(nrow(smokj) > 50)
	{
		gwas_list <- c(gwas_list, paste0(basename(out_file), "_lt25"))
		m <- "GWAS will be performed"
	} else {
		m <- "GWAS won't be performed"
	}
	message(nrow(smokj), " individuals under age 25 with predicted smoking data, ", m)
	write.table(smokj, file=paste0(out_file, "_lt25.plink"), row=F, col=F, qu=F)

	write.table(gwas_list, file=gwas_list_file, row=F, col=F, qu=F)
}



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

	dat <- data.frame(IID = colnames(mbeta), Smoking = scores_combined)
	return(dat)
}

rntransform <- function(x)
{
	out <- rank(x) - 0.5
	out[is.na(x)] <- NA
	mP <- 0.5/max(out, na.rm = T)
	out <- out/(max(out, na.rm = T) + 0.5)
	out <- scale(qnorm(out))
	out
}


main()