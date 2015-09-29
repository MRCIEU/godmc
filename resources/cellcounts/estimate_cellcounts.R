source('resources/cellcounts/wbcInference-V112.R')


##

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]

message("Reading methylation data...")
load(methylationfile)

cellcounts <- estimate.cellcounts(norm.beta)
cellcounts <- data.frame(IID=rownames(cellcounts), cellcounts)

write.table(cellcounts, file=cellcountfile, row=F, col=T, qu=F)