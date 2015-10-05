source('resources/cellcounts/wbcInference-V112.R')


##

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
ccrnmethdatafile <- arguments[3]
nthreads <- as.numeric(arguments[4])
chunks <- as.numeric(arguments[5])
jid <- as.numeric(arguments[6])


message("Reading methylation data...")
load(methylationfile)

cellcounts <- read.table(cellcountfile, he=T)
stopifnot(all(cellcounts$IID == colnames(norm.beta)))
cellcounts <- as.matrix(subset(cellcounts, select=-c(IID)))

if(!is.na(jid))
{
	i1 <- chunks * (jid-1) + 1
	i2 <- min(nrow(norm.beta), chunk * jid)
	norm.beta <- norm.beta[i1:i2,]
	ccrnmethdatafile <- paste0(ccrnmethdatafile, ".", jid, ".RData")
} else {
	ccrnmethdatafile <- paste0(ccrnmethdatafile, ".RData")
}

# adjust for cell counts
norm.beta <- adjust.beta(norm.beta, cellcounts, mc.cores=nthreads)

# and rank transformed data of cell countadjusted betas
norm.beta <- inverse.rank.transform(norm.beta)

# Save results
save(norm.beta, file=ccrnmethdatafile)
