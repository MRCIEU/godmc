source('resources/cellcounts/wbcInference-V112.R')


##

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
pcfile <- arguments[3]
ccrnmethdatafile <- arguments[4]
nthreads <- as.numeric(arguments[5])
chunks <- as.numeric(arguments[6])
jid <- as.numeric(arguments[7])


message("Reading methylation data...")
load(methylationfile)

cellcounts <- read.table(cellcountfile, he=T)
stopifnot(all(cellcounts$IID == colnames(norm.beta)))
cellcounts <- as.matrix(subset(cellcounts, select=-c(IID)))

load(pcfile)

if(!is.na(jid))
{
	chunksize <- ceiling(nrow(norm.beta) / chunks)
	i1 <- chunksize * (jid-1) + 1
	i2 <- min(nrow(norm.beta), chunksize * jid)
	norm.beta <- norm.beta[i1:i2,]
	ccrnmethdatafile <- paste0(ccrnmethdatafile, ".", jid, ".RData")
} else {
	ccrnmethdatafile <- paste0(ccrnmethdatafile, ".RData")
}

message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

# adjust for cell counts
norm.beta <- adjust.beta(norm.beta, pc, cellcounts, mc.cores=nthreads)

# and rank transformed data of cell countadjusted betas
norm.beta <- inverse.rank.transform(norm.beta, nthreads)

# Save results
save(norm.beta, file=ccrnmethdatafile)
