library(data.table)

#Thanks to EasyQC harmonization

args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
SNPfail.out <- as.character(args[2]);

#bim_file="./processed_data/genetic_data/data.bim"

bim <- as.data.frame(fread(bim_file))

#bim<-data.frame(f$V4,f[,1],f$V2,f$V3,f$V2,f$V3)
bim[,2]<-as.character(bim[,2])
bim[,5]<-as.character(bim[,5])
bim[,6]<-as.character(bim[,6])


harmonize.alleles <- function(bim=bim,SNPfail=SNPfail) {

message("Checking allele coding")
a1<-bim[,5]
a2<-bim[,6]

SNPfail<-data.frame()

### A1=NA & A2=NA
isBothMissing <- which(is.na(a1) & is.na(a2))
if(length(isBothMissing)>0) {SNPfail<-rbind(SNPfail,bim[isBothMissing,2])}

message("Allele harmonization:",length(isBothMissing)," alleles with  NA are going to be removed")
	
## Recode missing single allele or '<DEL>' or '-' to 'D' and set the other to 'I'

    a1<-bim[,5]
	a2<-bim[,6]
    isRecodea1 <- is.na(a1)|a1=='<DEL>'|a1=='-'
	bim[isRecodea1,5] <- "D"
	bim[isRecodea1,6] <- "I"
	
	message("Allele harmonization:",length(which(isRecodea1))," A1 alleles with  NA or <DEL> or - are recoded to D/I")
    a1<-bim[,5]
	a2<-bim[,6]
	isRecodea2 <- is.na(a2)|a2=='<DEL>'|a2=='-'
	bim[isRecodea2,5] <- "I"
	bim[isRecodea2,6] <- "D"
	
	message("Allele harmonization:",length(which(isRecodea2))," A2 alleles with  NA or <DEL> or - are recoded to I/D")
	
## Recode MACH R/I -> D/I and R/D -> I/D
	a1<-bim[,5]
	a2<-bim[,6]
    isRecodea1_mach1 <- a1=="R" & a2=="I"
	bim[isRecodea1_mach1,5] <- "D"
	message("Allele harmonization:",length(which(isRecodea1_mach1))," A1 alleles with  R/I are recoded to D/I")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecodea1_mach2 <- a1=="R" & a2=="D"
	bim[isRecodea1_mach2,5] <- "I"
	message("Allele harmonization:",length(which(isRecodea1_mach2))," A1 alleles with  R/D are recoded to I/D")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecodea2_mach1 <- a2=="R" & a1=="I"
	bim[isRecodea2_mach1,6] <- "D"
	message("Allele harmonization:",length(which(isRecodea2_mach1))," A2 alleles with  I/R are recoded to I/D")
	
	a1<-bim[,5]
    a2<-bim[,6]
	isRecodea2_mach2 <- a2=="R" & a1=="D"
	bim[isRecodea2_mach2,6] <- "I"
	message("Allele harmonization:",length(which(isRecodea2_mach2))," A2 alleles with  D/R are recoded to D/I")
	
## Recode Sequence coding to D/I
	a1<-bim[,5]
	a2<-bim[,6]
	isRecode_seq1 <- nchar(a1)>nchar(a2)
	bim[isRecode_seq1,5] <- "I"
	bim[isRecode_seq1,6] <- "D"
	message("Allele harmonization:",length(which(isRecode_seq1))," alleles with sequence coding are recoded to I/D")
	
	a1<-bim[,5]
	a2<-bim[,6]
	isRecode_seq2 <- nchar(a1)<nchar(a2)
	bim[isRecode_seq2,5] <- "D"
	bim[isRecode_seq2,6] <- "I"
	message("Allele harmonization:",length(which(isRecode_seq2))," alleles with sequence coding are recoded to D/I")

#Some SNPs in phase3 are coded as sequences
	#a1<-bim[,5]
	#a2<-bim[,6]
	#isRecode_seq3<- nchar(a1)>1 & nchar(a2)>1&nchar(a1)==nchar(a2)
	#bim[isRecode_seq3,2]<-gsub("INDEL","SNP",bim[isRecode_seq3,2])
    

   #testa1<-strsplit(as.character(bim[isRecode_seq3,5]),split="",fixed=T)
   #testa1<-unlist(lapply(testa1,function(x){y<-x[1]}))
   #testa2<-strsplit(as.character(bim[isRecode_seq3,6]),split="",fixed=T)
   #testa2<-unlist(lapply(testa2,function(x){y<-x[1]}))
   #w1<-testa1!=testa2
   #w2<-testa1==testa2
   #bim[which(isRecode_seq3)[w1],5] <-testa1[w1]
   #bim[which(isRecode_seq3)[w1],6] <-testa2[w1]
   #message("Allele harmonization:",length(which(isRecode_seq3)[w1])," SNPs with sequence coding are recoded to biallelic SNPs, for example AAGTTA/TAGTTA is recoded to A/T")
   
   #testa1<-strsplit(as.character(bim[which(isRecode_seq3)[w2],5]),split="",fixed=T)
   #testa1<-unlist(lapply(testa1,function(x){y<-x[length(x)]}))
   #testa2<-strsplit(as.character(bim[which(isRecode_seq3)[w2],6]),split="",fixed=T)
   #testa2<-unlist(lapply(testa2,function(x){y<-x[length(x)]}))
   #w3<-testa1!=testa2
   
   #bim[which(isRecode_seq3)[w2][w3],5] <-testa1[w3]
   #bim[which(isRecode_seq3)[w2][w3],6] <-testa2[w3]
   #message("Allele harmonization:",length(which(isRecode_seq3)[w2][w3])," SNPs with sequence coding are recoded to biallelic SNPs, for example AAATT/AAATA is recoded to T/A")
   
	a1<-bim[,5]
    a2<-bim[,6]
    isInvalid <- !(a1%in%c("A","C","G","T","I","D")&a2%in%c("A","C","G","T","I","D")&a1!=a2)
	if(length(which(isInvalid))>0) {
	SNPfail<-rbind(SNPfail,bim[which(isInvalid),2])}
    message("Allele harmonization:",length(which(isInvalid))," alleles with coding other than A,C,T,G,I,D are going to be removed")
	

	rm(a1,a2)
	rm(isRecode_seq1,isRecode_seq2)
	rm(isRecodea1_mach1,isRecodea1_mach2,isRecodea2_mach1,isRecodea2_mach2)
	rm(isInvalid)
	
	recoded.bim <- list(bim, SNPfail)
	
	return(recoded.bim)
}

recoded.bim<-harmonize.alleles(bim,SNPfail)
if(length(recoded.bim[[2]])>0){
SNPfailures<-as.character(t(recoded.bim[[2]]))
write.table(SNPfailures,SNPfail.out,sep="\t",quote=F,col.names=F,row.names=F)}
write.table(recoded.bim[[1]],bim_file,sep="\t",quote=F,col.names=F,row.names=F)

		
