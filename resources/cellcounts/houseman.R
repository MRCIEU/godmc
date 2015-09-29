# Adapted from https://gist.github.com/brentp/5058805

# Make a rank transformed methylation set
# Make a cell-count-adjusted and rank transformed methylation set
# These values will not be 0-1


source('resources/cellcounts/wbcInference-V112.R')


##

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
rnmethdatafile <- arguments[3]
ccrnmethdatafile <- arguments[4]
ccrnsquaredmethdatafile <- arguments[5]
nthreads <- as.numeric(arguments[6])

message("Reading methylation data...")
load(methylationfile)

if(cellcounts != "NULL")
{
	cellcounts <- read.table(cellcountfile, he=T)
	cellcounts <- subset(cellcounts, select=-c(IID))
}

# Get inverse rank transformed data, no cell count adjustment
dat <- inverse.rank.transform(norm.beta, nthreads)
write.table(round(dat, 3), file=rnmethdatafile, row=TRUE, col=TRUE, qu=FALSE, sep="\t")

# adjust for cell counts
if(cellcounts != "NULL")
{
	dat <- adjust.beta(norm.beta, cellcounts, mc.cores=nthreads)
	# and rank transformed data of cell countadjusted betas
	dat <- inverse.rank.transform(dat)
	write.table(round(dat, 3), file=ccrnmethdatafile, row=TRUE, col=TRUE, qu=FALSE, sep="\t")
} else {
	system(paste("cp -v", rnmethdatafile, ccrnmethdatafile))
}

# Get squared z values of cell count adjusted data
dat <- dat^2
write.table(round(dat, 3), file=ccrnsquaredmethdatafile, row=TRUE, col=TRUE, qu=FALSE, sep="\t")
