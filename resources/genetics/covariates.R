# Create covariates file

# Original covariates
# Genetic PCs
# Cell counts


arguments <- commandArgs(T)

cov_file <- arguments[1]
pca_file <- arguments[2]
cellcount_file <- arguments[3]
out_file <- arguments[4]

pca <- read.table(pca_file)[,-1]
names(pca) <- c("IID", paste("genetic_pc", 1:(ncol(pca)-1), sep=""))
if(cov_file=="NULL")
{
	covs <- data.frame(IID=pca$IID)
} else {
	covs <- read.table(cov_file, header=TRUE)
}

cellcount <- read.table(cellcount_file, header=TRUE)

allcovs <- merge(cellcount, pca, by="IID", all=TRUE)
allcovs <- merge(allcovs, covs, by="IID", all=TRUE)

mat <- t(as.matrix(allcovs))
write.table(mat, file=out_file, row=TRUE, col=TRUE, qu=FALSE, sep="\t")

