arguments <- commandArgs(T)

betas <- arguments[1]
cellcountfile <- arguments[2]
cellcountref <- arguments[3]
#cellcountplinkraw <- arguments[4]
fam_file<- arguments[4]

cellcountref<-gsub("_"," ",cellcountref)
cat(cellcountref,"\n")

library(meffil)
load(betas)
cellcounts<-meffil.estimate.cell.counts.from.betas(norm.beta,cell.type.reference=cellcountref,verbose=T)
cellcounts<-data.frame(IID=row.names(cellcounts),cellcounts)
fam <- read.table(fam_file, stringsAsFactors=FALSE)
m<-match(fam[,2],cellcounts$IID)
write.table(cellcounts[m,],cellcountfile,sep="\t",row.names=F,col.names=T,quote=F)
#cellcounts<-data.frame(FID=fam[,1],cellcounts[m,])
#write.table(cellcounts,cellcountplinkraw,sep="\t",row.names=F,col.names=F,quote=F)
