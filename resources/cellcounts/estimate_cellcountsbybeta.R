suppressMessages(library(meffil))


arguments <- commandArgs(T)

betas <- arguments[1]
cellcountfile <- arguments[2]
cellcountref <- arguments[3]
fam_file <- arguments[4]


message("Loading methylation betas")
load(betas)

cellcountref <- gsub("_"," ",cellcountref)
message("Using '", cellcountref, "' reference")

cellcounts <- meffil.estimate.cell.counts.from.betas(norm.beta,cell.type.reference=cellcountref,verbose=T)
cellcounts <- data.frame(IID=row.names(cellcounts),cellcounts)
fam <- read.table(fam_file, stringsAsFactors=FALSE)
m <- match(fam[,2],cellcounts$IID)
write.table(cellcounts[m,],cellcountfile,sep="\t",row.names=F,col.names=T,quote=F)
