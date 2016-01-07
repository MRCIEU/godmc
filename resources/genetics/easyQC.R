args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
frq_file <- as.character(args[2]);
out_file <- as.character(args[3]);
easyQC_file <- as.character(args[4]);
easyQC_script <- as.character(args[5]);

library(data.table)
bim<-as.data.frame(fread(bim_file))
message("read bimfile")
frq<-as.data.frame(fread(frq_file,header=T))
message("read frq file")

m<-match(bim[,2],frq[,2])
bim<-data.frame(bim,frq$MAF,N=frq$NCHROBS/2,BETA=0,SE=0,PVAL=1,IMPUTATION=1,STRAND="+")
names(bim)<-c("CHR","SNP","GeneticPOS","POS","EFFECT_ALLELE","OTHER_ALLELE","EAF","N","BETA","SE","PVAL","IMPUTATION","STRAND")
write.table(bim[,-3],paste(out_file,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
message("created easyQC input file")

library(EasyQC)
EasyQC(easyQC_script)
message("easyQC strand check finished successfully")
message("read easyQC_file")
r<-read.table(easyQC_file,header=T,sep="\t")
cat(easyQC_file,"\n")
head(easyQC_file)
message("read easyQC file successfully")

r<-r[match(bim[,2],r$SNP),]
message("read mismatch file")
mm<-read.table(paste(out_file,".mismatch.txt",sep=""),header=T)
message("read mismatch file successfully")

message("read not in ref file")
notinref<-read.table(paste(out_file,".notinref.txt",sep=""),header=T)
message("read not in ref file successfully")

message("read AF outlier file")
afoutlier<-read.table(paste(out_file,".AFCHECK.outlier.txt",sep=""),sep="\t",header=T)
message("read AF outlier file successfully")

strand.rm<-unique(c(as.character(mm$SNP),as.character(notinref$SNP),as.character(afoutlier$SNP)))
w<-which(bim[,2]%in%strand.rm) #4494

a1<-which(bim[,5]!=r[,3] & bim[,6]!=r[,4]) #106629
a1<-a1[a1%in%w==F] #102135
#4494

SNPflip<-as.character(bim[a1,2])
write.table(SNPflip,paste(out_file,".flipped.SNPs.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)

SNP.rm<-unique(c(as.character(mm$SNP),as.character(afoutlier$SNP)))
write.table(SNP.rm,paste(out_file,".mismatch_afcheck.failed.SNPs.txt",sep=""),sep="\t",quote=F,row.names=F,col.names=F)
message("extracted mismatched and misaligned SNPs")