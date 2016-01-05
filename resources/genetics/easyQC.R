args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
frq_file <- as.character(args[2]);
out_file <- as.character(args[3]);
easyQC_file <- as.character(args[4]);
easyQC_script <- as.character(args[5]);

library(data.table)

bim<-read.table(paste("./",bim_file,sep=""))
message("read bimfile")
frq<-read.table(paste("./",frq_file,sep=""),header=T)
message("read frq file")

m<-match(bim[,2],frq[,2])
bim<-data.frame(bim,frq$MAF,N=frq$NCHROBS/2,BETA=0,SE=0,PVAL=1,IMPUTATION=1,STRAND="+")
names(bim)<-c("CHR","SNP","GeneticPOS","POS","EFFECT_ALLELE","OTHER_ALLELE","EAF","N","BETA","SE","PVAL","IMPUTATION","STRAND")
write.table(bim[,-3],paste("./",out_file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
message("created easyQC input file")


library(EasyQC)
EasyQC(paste("./",easyQC_script,sep="\t")
message("easyQC strand check finished successfully")

r<-read.table(paste("./",easyQC_file,sep=""),header=T)
message("read easyQC file")

r<-r[match(bim[,2],r$SNP),]
mm<-read.table(paste("./",out_file,".mismatch.txt",sep=""),header=T)
notinref<-read.table(paste("./",out_file,".notinref.txt",sep=""),header=T)
afoutlier<-read.table(paste("./",out_file,".AFCHECK.outlier.txt",sep=""),header=T)

strand.rm<-unique(c(as.character(mm$SNP),as.character(notinref$SNP),as.character(afoutlier$SNP)))
w<-which(bim[,2]%in%strand.rm) #4494

a1<-which(bim[,5]!=r[,3] & bim[,6]!=r[,4]) #106629
a1<-a1[a1%in%w==F] #102135
#4494

SNPflip<-as.character(bim[a1,2])
write.table(SNPflip,paste("./",out_file,".flipped.SNPs.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

SNP.rm<-unique(c(as.character(mm$SNP),as.character(afoutlier$SNP)))
write.table(SNP.rm,paste("./",out_file,".mismatch_afcheck.failed.SNPs.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
message("extracted mismatched and misaligned SNPs")