source('resources/cellcounts/wbcInference-V112.R')


##

arguments <- commandArgs(T)

methylationfile <- arguments[1]
rnmethdatafile <- arguments[2]
nthreads <- as.numeric(arguments[3])

##

message("Reading methylation data...")
load(methylationfile)

# Get inverse rank transformed data
norm.beta <- inverse.rank.transform(norm.beta, nthreads)
save(norm.beta, file=paste0(rnmethdatafile, ".RData"))
