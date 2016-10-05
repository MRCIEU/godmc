	arguments <- commandArgs(T)

	gctafile <- arguments[1]
	plinkfile <- arguments[2]
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

pdf(paste(outdir,"/plinkvsgcta.",probeset,".pdf",sep=""),height=6,width=6)
plot(-log10(plink$P),-log10(gcta$P),cex=0.8,xlab="-log10 P (plink)",ylab="-log10 P (gcta)")
dev.off()

plinkP<--log10(plink$P)
gctaP<--log10(gcta$P)
diff<-as.numeric(gctaP)-as.numeric(plinkP)
diff.abs<-abs(diff)
head(plink[which(diff.abs>50),c("id","EA","EAF","BETA","P")])
head(gcta[which(diff.abs>50),c("id","EA","EAF","BETA","P")])