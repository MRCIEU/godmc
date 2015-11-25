    args = (commandArgs(TRUE));
    cellcount_file = as.character(args[1]);
	out_file = as.character(args[2]);
	ids = as.character(args[3]);
    cellcounts.plot = as.character(args[4]);
    covfile = as.character(args[5]);
    cellcounts.plink = as.character(args[6]);
    SD = as.numeric(args[7]);
    transformed.cellcounts = as.character(args[8]);
    cellcounts.entropy = as.character(args[9]);
    smoking.pred = as.character(args[10]);
    transformed.cellcounts.smokadj = as.character(args[11]);
    cellcounts.smokadj.plink = as.character(args[12]);
#cellcount_file <- "./input_data/cellcounts.txt"
#ids<-("./processed_data/ids/ids_plink.txt")
#cellcounts.plot<-"cellcountsplots.pdf"
#covfile<-"./input_data/covariates.txt"
#transformed.cellcounts<-"transformed.cellcounts.txt"
#SD=5

library(lattice)

#
#fam_file<-"GoDMC_data.fam"
#out_file<-"rawcellcountsforGWAs.txt"
	
	cellcounts <- read.table(cellcount_file, he=T)
	#fam <- read.table(fam_file)[,1:2]

	cellcounts$entropy <- apply(as.matrix(cellcounts[,-1]), 1, function(x)
	{
		x <- x[x>0]
		h <- x * log2(x)
		h[!is.finite(h)] <- 0
		-sum(h)
	})

    IID<-read.table(ids)
    m<-match(IID[,2],cellcounts$IID)
    cellcounts<-cellcounts[m,]
    write.table(cellcounts[m,],file=cellcounts.entropy,sep="\t",row.names=F,col.names=T,quote=F)
    cellcounts.out<-data.frame(FID=IID[,1],cellcounts[m,])
    write.table(cellcounts.out, file=out_file, row=F, col=F, qu=F)

#	for(i in 2:ncol(cellcounts))
#	{
#		cellcounts[,i] <- rntransform(cellcounts[,i])
#	}
#	nom <- names(cellcounts)[-1]
#	cellcounts <- merge(cellcounts, fam, by.x="IID", by.y="V2")
#	cellcounts <- subset(cellcounts, select=c(V1, V2, nom))

#	write.table(cellcounts, file=out_file, row=F, col=F, qu=F)	
#}


traits<-names(cellcounts)[-1]
outdata.all<-as.character(IID[,2])
outdata.all2<-as.character(IID[,2])



pdf(cellcounts.plot, width=12, height=8)
par(mfrow=c(2,2))

for (tr in 1:length(traits)){
data<- cellcounts
cov<-read.table(covfile,header=T,stringsAsFactors =F)
m<-match(data$IID,cov$IID)
data<-data.frame(data,cov[m,-1])
smoking<-read.table(paste(smoking.pred,".txt",sep=""),header=T)
m<-match(data$IID,smoking$IID)
data<-data.frame(data,cov[m,-1],Smoking=smoking[m,-1])

trait_var =traits[tr]

data$trait <- data[[trait_var]]
#data <- data[which(data$Sex!="NA" &data$trait!="NA"& data$trait!="na"&data$trait!="0"),] 


data$trait <- as.numeric(data$trait)

w.test=1
if (length(which(names(data)%in%"Sex"))>0 & length(rownames(table(data$Sex)))>1){
data <- data[which(data$Sex!="NA" &data$trait!="NA"& data$trait!="na"),] 
data$Sex <- as.factor(data$Sex)

sex_names <- rownames(table(data$Sex))
w.test<-wilcox.test(data$trait[data$Sex=="M"], data$trait[data$Sex=="F"] )[3]

colors <- 1:length(sex_names)
plot(density( data$trait[which(data$Sex==sex_names[1])],na.rm=T),xlab="",main=paste(trait_var," density plot by Sex",sep=""), col=colors[1])
for (j in 2:3) { ## a maximum of 3 types of Sex including unknown
	subset <- data$trait[which(data$Sex==sex_names[j] & !is.na(data$trait))]
	if (length(sex_names) >=j & length(subset) >0) {
		lines(density(subset,na.rm=T),xlab="",col=colors[j], main="")
	}
}
legend("topright",sex_names,col=colors,lty=1)	
}

#### plot the distribution of raw phenotypes
par(mfrow=c(2,2))
if (w.test>0.05|length(rownames(table(data$Sex)))==1){
data <- data[which(data$trait!="NA"& data$trait!="na"),]
plot(data$trait, xlab="", main=paste("raw ",trait_var," (N=", length(which(!is.na(data$trait))),")",sep=""),cex.main=0.7)
hist(data$trait, xlab="", main=paste("raw ",trait_var," (N=", length(which(!is.na(data$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(data$trait,na.rm=T)-SD*sd(data$trait,na.rm=T),lty=2)
abline(v=mean(data$trait,na.rm=T)+SD*sd(data$trait,na.rm=T),lty=2)
qqnorm(data$trait, main=paste("raw ",trait_var," (N=", length(which(!is.na(data$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(data$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(data$trait)

par(mfrow=c(2,2))

#remove outlier
outlier <- which(data$trait<(mean(data$trait,na.rm=T)-SD*sd(data$trait,na.rm=T)) | data$trait> (mean(data$trait,na.rm=T)+SD*sd(data$trait,na.rm=T)))
if (length(outlier)>0){data<-data[-outlier,]}

#transform data
data$trait <- qnorm((rank(data$trait,na.last="keep")-0.5)/sum(!is.na(data$trait)))
data$trait_smokadj<-data$trait
#adjust for age
if(length(which(names(data)%in%c("Age")))==1){

fit1<- lm(trait ~ Age, data=data)
fit2<- lm(trait ~ Age+I(Age^2), data=data)
fit3<- lm(trait ~ Age + Smoking, data=data)
fit4<- lm(trait ~ Age+I(Age^2)+Smoking, data=data)
fit5<- lm(trait ~ Smoking, data=data)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]<0.05){
data$trait<-resid(fit1)
data$trait_smokadj<-resid(fit3)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]<0.05){
data$trait<-resid(fit2)
data$trait_smokadj<-resid(fit4)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]>0.05 & coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]>0.05){
data$trait_smokadj<-resid(fit5)
}

if(length(which(names(data)%in%c("Age")))!=1){
data$trait_smokadj<-resid(fit5)
}

}
#standardise
nmiss<-which(!is.na(data[,"trait"]))
data[nmiss,"trait"]<-(data[nmiss,"trait"]-mean(data[nmiss,"trait"]))/sd(data[nmiss,"trait"])

nmiss<-which(!is.na(data[,"trait_smokadj"]))
data[nmiss,"trait_smokadj"]<-(data[nmiss,"trait_smokadj"]-mean(data[nmiss,"trait_smokadj"]))/sd(data[nmiss,"trait_smokadj"])

outdata<-data.frame(IID=data$IID,trait=data$trait,trait_smokadj=data$trait_smokadj)
m<-match(IID[,2],outdata$IID)
outdata<-data.frame(IID=IID[,2],trait=outdata$trait[m],trait_smokadj=outdata$trait_smokadj[m])

#### plot the distribution of transformed phenotypes
par(mfrow=c(2,2))
plot(outdata$trait, xlab="", main=paste("transformed ",trait_var," (N=", length(which(!is.na(outdata$trait))),")",sep=""),cex.main=0.7)
hist(outdata$trait, xlab="", main=paste("transformed ",trait_var," (N=", length(which(!is.na(outdata$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(outdata$trait,na.rm=T)-SD*sd(outdata$trait,na.rm=T),lty=2)
abline(v=mean(outdata$trait,na.rm=T)+SD*sd(outdata$trait,na.rm=T),lty=2)
qqnorm(outdata$trait, main=paste("transformed ",trait_var," (N=", length(which(!is.na(outdata$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(outdata$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(outdata$trait)
par(mfrow=c(2,2))

}else{

male <- data[which(data$Sex=="M"),]
female <- data[which(data$Sex=="F"),]

par(mfrow=c(2,2))
plot(male$trait, xlab="", main=paste("raw ",trait_var," (males, N=", length(which(!is.na(male$trait))),")",sep=""),cex.main=0.7)
hist(male$trait, xlab="", main=paste("raw ",trait_var," (males, N=", length(which(!is.na(male$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(male$trait,na.rm=T)-SD*sd(male$trait,na.rm=T),lty=2)
abline(v=mean(male$trait,na.rm=T)+SD*sd(male$trait,na.rm=T),lty=2)
qqnorm(male$trait, main=paste("raw ",trait_var," (males, N=", length(which(!is.na(male$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(male$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(male$trait)

par(mfrow=c(2,2))
plot(female$trait, xlab="", main=paste("raw ",trait_var," (females, N=", length(which(!is.na(female$trait))),")",sep=""),cex.main=0.7)
hist(female$trait, xlab="", main=paste("raw ",trait_var," (females, N=", length(which(!is.na(female$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(female$trait,na.rm=T)-SD*sd(female$trait,na.rm=T),lty=2)
abline(v=mean(female$trait,na.rm=T)+SD*sd(female$trait,na.rm=T),lty=2)
qqnorm(female$trait, main=paste("raw ",trait_var," (females, N=", length(which(!is.na(female$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(female$trait)[2]),2),")",sep=""),cex.main=0.6)
qqline(female$trait)

#remove outliers
outlierm <- which(male$trait<(mean(male$trait,na.rm=T)-SD*sd(male$trait,na.rm=T)) | male$trait> (mean(male$trait,na.rm=T)+SD*sd(male$trait,na.rm=T)))
if (length(outlierm)>0){male<-male[-outlierm,]}
outlierf <- which(female$trait<(mean(female$trait,na.rm=T)-SD*sd(female$trait,na.rm=T)) | female$trait> (mean(female$trait,na.rm=T)+SD*sd(female$trait,na.rm=T)))
if (length(outlierf)>0){female<-female[-outlierf,]}

#transform
male$trait <-qnorm((rank(male$trait,na.last="keep")-0.5)/sum(!is.na(male$trait)))
female$trait<-qnorm((rank(female$trait,na.last="keep")-0.5)/sum(!is.na(female$trait)))
male$trait_smokadj<-male$trait
female$trait_smokadj<-female$trait
#adjust for covariates

if(length(which(names(data)%in%"Age"))>0){

fit1<- lm(trait ~ Age, data=male)
fit2<- lm(trait ~ Age+I(Age^2), data=male)
fit3<- lm(trait ~ Age + Smoking, data=male)
fit4<- lm(trait ~ Age+I(Age^2)+Smoking, data=male)
fit5<- lm(trait ~ Smoking, data=male)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]<0.05){
male$trait<-resid(fit1)
male$trait_smokadj<-resid(fit3)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]<0.05){
male$trait<-resid(fit2)
male$trait_smokadj<-resid(fit4)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]>0.05 & coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]>0.05){
male$trait_smokadj<-resid(fit5)
}

if(length(which(names(male)%in%c("Age")))!=1){
male$trait_smokadj<-resid(fit5)
}

fit1<- lm(trait ~ Age, data=female)
fit2<- lm(trait ~ Age+I(Age^2), data=female)
fit3<- lm(trait ~ Age + Smoking, data=female)
fit4<- lm(trait ~ Age+I(Age^2)+Smoking, data=female)
fit5<- lm(trait ~ Smoking, data=female)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]<0.05){
female$trait<-resid(fit1)
female$trait_smokadj<-resid(fit3)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]<0.05){
female$trait<-resid(fit2)
female$trait_smokadj<-resid(fit4)
}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age^2)"]>0.05 & coefficients(summary(fit1))[,"Pr(>|t|)"]["Age"]>0.05){
female$trait_smokadj<-resid(fit5)
}

if(length(which(names(female)%in%c("Age")))!=1){
female$trait_smokadj<-resid(fit5)
}
}

#standardise
nmiss_male<-which(!is.na(male[,"trait"]))
male[nmiss_male,"trait"]<-(male[nmiss_male,"trait"]-mean(male[nmiss_male,"trait"]))/sd(male[nmiss_male,"trait"])
nmiss_female<-which(!is.na(female[,"trait"]))
female[nmiss_female,"trait"]<-(female[nmiss_female,"trait"]-mean(female[nmiss_female,"trait"]))/sd(female[nmiss_female,"trait"])


nmiss_male<-which(!is.na(male[,"trait_smokadj"]))
male[nmiss_male,"trait_smokadj"]<-(male[nmiss_male,"trait_smokadj"]-mean(male[nmiss_male,"trait_smokadj"]))/sd(male[nmiss_male,"trait_smokadj"])
nmiss_female<-which(!is.na(female[,"trait_smokadj"]))
female[nmiss_female,"trait_smokadj"]<-(female[nmiss_female,"trait_smokadj"]-mean(female[nmiss_female,"trait_smokadj"]))/sd(female[nmiss_female,"trait_smokadj"])

outdata<-rbind(male,female)
m<-match(IID[,2],outdata$IID)
outdata<-data.frame(IID=IID[,2],trait=outdata$trait[m],trait_smokadj=outdata$trait_smokadj[m])

#after transformation
par(mfrow=c(2,2))
plot(male$trait, xlab="", main=paste("transformed ",trait_var," (males, N=", length(which(!is.na(male$trait))),")",sep=""),cex.main=0.7)
hist(male$trait, xlab="", main=paste("transformed ",trait_var," (males, N=", length(which(!is.na(male$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(male$trait,na.rm=T)-SD*sd(male$trait,na.rm=T),lty=2)
abline(v=mean(male$trait,na.rm=T)+SD*sd(male$trait,na.rm=T),lty=2)
qqnorm(male$trait, main=paste("transformed ",trait_var," (males, N=", length(which(!is.na(male$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(male$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(male$trait)

par(mfrow=c(2,2))
plot(female$trait, xlab="", main=paste("transformed ",trait_var," (females, N=", length(which(!is.na(female$trait))),")",sep=""),cex.main=0.7)
hist(female$trait, xlab="", main=paste("transformed ",trait_var," (females, N=", length(which(!is.na(female$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(female$trait,na.rm=T)-SD*sd(female$trait,na.rm=T),lty=2)
abline(v=mean(female$trait,na.rm=T)+SD*sd(female$trait,na.rm=T),lty=2)
qqnorm(female$trait, main=paste("transformed ",trait_var," (females, N=", length(which(!is.na(female$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(female$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(female$trait)

par(mfrow=c(2,2))
plot(outdata$trait, xlab="", main=paste("transformed combined ",trait_var," (N=", length(which(!is.na(outdata$trait))),")",sep=""),cex.main=0.7)
hist(outdata$trait, xlab="", main=paste("transformed combined ",trait_var," (N=", length(which(!is.na(outdata$trait))),")",sep=""),cex.main=0.7)
abline(v=mean(outdata$trait,na.rm=T)-SD*sd(outdata$trait,na.rm=T),lty=2)
abline(v=mean(outdata$trait,na.rm=T)+SD*sd(outdata$trait,na.rm=T),lty=2)
qqnorm(outdata$trait, main=paste("transformed combined ",trait_var," (N=", length(which(!is.na(outdata$trait))),"; shapiroP=",signif(as.numeric(shapiro.test(outdata$trait)[2]),2),")",sep=""),cex.main=0.7)
qqline(outdata$trait)
par(mfrow=c(2,2))

}
outdata.all<-cbind(outdata.all,outdata$trait)
outdata.all2<-cbind(outdata.all2,outdata$trait_smokadj)
}
null <- dev.off()

colnames(outdata.all)<-c("IID",traits)
colnames(outdata.all2)<-c("IID",traits)
write.table(outdata.all,transformed.cellcounts,sep="\t",quote=F,row.names=F,col.names=T)
write.table(outdata.all2,transformed.cellcounts.smokadj,sep="\t",quote=F,row.names=F,col.names=T)


m<-match(IID[,2],outdata.all[,1])
outdata.all<-data.frame(FID=IID[,1],outdata.all[m,])
# names(outdata.all)

write.table(outdata.all,file=cellcounts.plink,sep="\t",quote=F,row.names=F,col.names=F)

# names(outdata.all2)
m<-match(IID[,2],outdata.all2[,1])
outdata.all2<-data.frame(FID=IID[,1],outdata.all2[m,])
# names(outdata.all2)

write.table(outdata.all2,file=cellcounts.smokadj.plink,sep="\t",quote=F,row.names=F,col.names=F)


