library(meffil)

arguments <- commandArgs(T)

methylationfile <- arguments[1]
cellcountfile <- arguments[2]
cellcountref <- arguments[3]
fam_file<- arguments[4]

message("Reading methylation data...")
load(methylationfile)

#index <- grepl("rs", rownames(norm.beta))
#cellcounts <- meffil.estimate.cell.counts.from.betas(norm.beta[!index, ], cellcountref)

cellcounts<-meffil.estimate.cell.counts.from.betas(norm.beta,cell.type.reference=cellcountref,verbose=T)
cellcounts<-data.frame(IID=row.names(cellcounts),cellcounts)
fam <- read.table(fam_file, stringsAsFactors=FALSE)
m<-match(fam[,2],cellcounts$IID)
write.table(cellcounts[m,],cellcountfile,sep="\t",row.names=F,col.names=T,quote=F)
