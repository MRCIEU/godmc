library(meffil)

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
cellcountref <- arguments[3]


message("Reading methylation data...")
load(methylationfile)

index <- grepl("rs", rownames(norm.beta))
cellcounts <- meffil.estimate.cell.counts.from.betas(norm.beta[!index, ], cellcountref)
cellcounts <- data.frame(IID=rownames(cellcounts), cellcounts)

write.table(cellcounts, file=cellcountfile, row=F, col=T, qu=F)