	arguments <- commandArgs(T)

	plinkfile <- arguments[1]
	gctafile <- arguments[2]
	probeset <- arguments[3]
    outdir<-arguments[4]
    
message("Reading plink data...")
plink<-read.table(plinkfile,he=T)
plink<-data.frame(id=paste(plink$CpG,plink$SNP,sep="_"),plink)

message("Reading gcta data...")
gcta<-read.table(gctafile,he=T)
gcta<-data.frame(id=paste(gcta$CpG,gcta$SNP,sep="_"),gcta)

nr1<-nrow(plink)
nr2<-nrow(gcta)

if(nr1!=nr2)
{
	msg <- paste0("ERROR:plink file has different number of rows  compared to gcta file")
	
}

m<-match(gcta$id,plink$id)
plink<-plink[m,]

diff<-(as.numeric(-log10(gcta$P)))-(as.numeric(-log10(plink$P)))
diff.abs<-abs(diff)
w<-which(diff.abs>10) #1571/18535=0.08475856
probe<-unique((plink[w,"CpG"]))
cols=rainbow(length(probe))

m<-match(plink$CpG,probe)
cols<-cols[m]

cols2<-cols
b1<-which(is.na(cols))
cols2[b1]<-"black"

pdf(paste(outdir,"/plinkvsgcta.",probeset,".pdf",sep=""),height=6,width=6)
plot(-log10(plink$P),-log10(gcta$P),cex=0.7,xlab="-log10 P (phase 1 - plink)",ylab="-log10 P (phase 2 - gcta)",pch=16,col=cols2,main=paste(length(unique(plink$CpG))," probes",sep=""))
abline(a=0,b=1,col="grey")
legend("bottomright",legend=probe,col=na.omit(unique(cols)),cex=0.5,lty=1, title="probes with -log10 P diff >10")

plot(-log10(plink$P),-log10(gcta$P),cex=0.7,xlab="-log10 P (phase 1 - plink)",ylab="-log10 P (phase 2 - gcta)",pch=16,col=cols2,main=paste(length(unique(plink$CpG))," probes",sep=""),ylim=c(0,13),xlim=c(0,13))
abline(a=0,b=1,col="grey")
legend("bottomright",legend=probe,col=na.omit(unique(cols)),cex=0.5,lty=1, title="probes with -log10 P diff >10")

plot(-log10(plink$P),-log10(gcta$P),cex=0.7,xlab="-log10 P (phase 1 - plink)",ylab="-log10 P (phase 2 - gcta)",pch=16,col=cols2,main=paste(length(unique(plink$CpG))," probes",sep=""),ylim=c(0,5),xlim=c(0,5))
abline(a=0,b=1,col="grey")
legend("bottomright",legend=probe,col=na.omit(unique(cols)),cex=0.5,lty=1, title="probes with -log10 P diff >10")


plot(plink$EAF,gcta$EAF,cex=0.7,xlab="Effect Allele Freq (plink)",ylab="Effect Allele Freq (gcta)",pch=16)
abline(a=0,b=1,col="grey")

plot(plink$BETA,gcta$BETA,cex=0.7,xlab="Effect Size (plink)",ylab="Effect Size (gcta)",pch=16,ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),col=cols2)
legend("bottomright",legend=probe,col=na.omit(unique(cols)),cex=0.5,lty=1, title="probes with -log10 P diff >10")
abline(a=0,b=1,col="grey")
dev.off()

