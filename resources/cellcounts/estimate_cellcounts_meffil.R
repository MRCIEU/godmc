arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
cellcountref <- arguments[3]


message("Reading methylation data...")
load(methylationfile)


cellcounts <- meffil.estimate.cellcounts.from.betas(norm.beta, cellcountref)
cellcounts <- data.frame(IID=rownames(cellcounts), cellcounts)

write.table(cellcounts, file=cellcountfile, row=F, col=T, qu=F)