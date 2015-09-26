args = (commandArgs(TRUE));
pathtoidat= as.character(args[1]);
cat(pathtoidat,"\n")

library("minfi")
#pathtoidat<-"/panfs/panasas01/sscm/epzjlm/ARIES/methylation/normalisation150508/test_idat"
Basename<-list.files(pathtoidat)
Basename<-gsub("_Grn.idat","",Basename)
Basename<-gsub("_Red.idat","",Basename)
Basename<-unique(Basename)

spl<-do.call("rbind",strsplit(Basename,"_"))
Basename<-paste(pathtoidat,Basename,sep="/")

targets<-data.frame(Basename=as.character(Basename))
rgset <- read.450k.exp(targets=targets)

class(rgset)
#[1] "RGChannelSet"

#make sure your columns are not coded as a factor.
rgset@phenoData@data$Basename<-as.character(rgset@phenoData@data$Basename)

#normalise your data together with cellcount data using preprocessQuantile...
pdf("cellcounts.pdf",height=6,width=6)
cellcounts<-estimateCellCounts(rgset,compositeCellType = "Blood",cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),returnAll = TRUE, meanPlot = TRUE, verbose = TRUE)
dev.off()

cellcounts<-data.frame(Sample_Name=row.names(cellcounts$counts),cellcounts$counts)
write.table(cellcounts,"cellcounts.txt",sep="\t",col.names=T,row.names=F,quote=F)

