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

diff<-as.numeric(gctaP)-as.numeric(plinkP)
diff.abs<-abs(diff)
w<-which(df$absPdiff>10) #1571/18535=0.08475856
probe<-unique((df[w,"CpG"]))
cols=rainbow(length(probe))
m<-match(df$CpG,probe)
cols<-cols[m]
b1<-which(is.na(cols))
cols[b1]<-"black"

pdf(paste(outdir,"/plinkvsgcta.",probeset,".pdf",sep=""),height=6,width=6)
plot(-log10(plink$P),-log10(gcta$P),cex=0.7,xlab="-log10 P (plink)",ylab="-log10 P (gcta)",pch=16,col=cols,main=paste(length(unique(plink$CpG))," probes",sep=""),)
abline(a=0,b=1,col="grey")
plot(-log10(plink$P),-log10(gcta$P),cex=0.7,xlab="-log10 P (plink)",ylab="-log10 P (gcta)",pch=16,col=cols,main=paste(length(unique(plink$CpG))," probes",sep=""),ylim=c(0,13),xlim=c(0,13))
abline(a=0,b=1,col="grey")

legend("bottomright",legend=probe,col=unique(cols),cex=0.5,lty=1)
plot(plink$EAF,gcta$EAF,cex=0.7,xlab="EAF (plink)",ylab="EAF (gcta)",pch=16)
abline(a=0,b=1,col="grey")
plot(plink$BETA,gcta$BETA,cex=0.7,xlab="Effect Size (plink)",ylab="Effect Size (gcta)",pch=16,ylim=c(-1.5,1.5),xlim=c(-1.5,1.5),col=cols)
legend("bottomright",legend=probe,col=unique(cols),cex=0.5,lty=1)
abline(a=0,b=1,col="grey")
dev.off()

plinkP<--log10(plink$P)
gctaP<--log10(gcta$P)
diff<-as.numeric(gctaP)-as.numeric(plinkP)
diff.abs<-abs(diff)
o<-order(diff.abs,decreasing=T)

head(plink[o,c("id","CpG","EA","EAF","BETA","P")])
head(gcta[o,c("id","CpG","EA","EAF","BETA","P")])