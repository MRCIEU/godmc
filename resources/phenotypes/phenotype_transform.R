args = (commandArgs(TRUE));
ids = as.character(args[1]);
phenofile= as.character(args[2]);
pheno.plot=as.character(args[3])
covfile= as.character(args[4]);
transformed.phenotypes=as.character(args[5]);
SD=as.numeric(args[6]);
phenotypes.plink=as.character(args[7]);
outfile=as.character(args[8]);

library(lattice)

#transformed.phenotypes<-"processed_data/phenotypes/transformed.phenotypes.txt"
#transformed.phenotypes<-"transformed.phenotypes.txt"
#phenofile <- "./input_data/EWAS.phenotypes.raw.txt"
#covfile<-"./input_data/covariates.txt"
#ids<-"./processed_data/ids/intersect_ids_plink.txt"
#covfile<-"covariates_withoutSex.txt"
#pheno.plot<-"phenotypeplots.pdf"
#SD=5

data<- read.table(phenofile, header=T,stringsAsFactors =F)
traits<-names(data)[-1]
IID<-read.table(ids)
outdata.all<-as.character(IID[,2])
m<-match(IID[,2],data$IID)
data2<-data.frame(FID=IID[,1],data[m,])
write.table(data2,file=paste(phenotypes.plink,".raw",sep=""),quote=F,row.names=F,col.names=F)
write.table(data2[,-1],file=paste0(outfile,".txt"),quote=F,row.names=F,col.names=T)

pdf(pheno.plot, width=12, height=8)
par(mfrow=c(2,2))

for (tr in 1:length(traits)){
data<- read.table(phenofile, header=T,stringsAsFactors =F)
cov<-read.table(covfile,header=T,stringsAsFactors =F)
m<-match(data$IID,cov$IID)
data<-data.frame(data,cov[m,-1])

trait_var =traits[tr]

data$trait <- data[[trait_var]]
#data <- data[which(data$Sex!="NA" &data$trait!="NA"& data$trait!="na"&data$trait!="0"),] 


data$trait <- as.numeric(data$trait)

w.test=1
if (length(which(names(data)%in%"Sex_factor"))>0 & length(rownames(table(data$Sex_factor)))>1){
data <- data[which(data$Sex_factor!="NA" &data$trait!="NA"& data$trait!="na"),] 
data$Sex_factor <- as.factor(data$Sex_factor)
sex_names <- rownames(table(data$Sex_factor))

w.test<-wilcox.test(data$trait[data$Sex_factor=="M"], data$trait[data$Sex_factor=="F"] )[3]

colors <- 1:length(sex_names)
plot(density( data$trait[which(data$Sex_factor==sex_names[1])],na.rm=T),xlab="",main=paste(trait_var," density plot by Sex",sep=""), col=colors[1])
for (j in 2:3) { ## a maximum of 3 types of Sex including unknown
	subset <- data$trait[which(data$Sex_factor==sex_names[j] & !is.na(data$trait))]
	if (length(sex_names) >=j & length(subset) >0) {
		lines(density(subset,na.rm=T),xlab="",col=colors[j], main="")
	}
}
legend("topright",sex_names,col=colors,lty=1)	
}

#### plot the distribution of raw phenotypes
par(mfrow=c(2,2))
if (w.test>0.05|length(rownames(table(data$Sex_factor)))==1){
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

#adjust for age
if(length(which(names(data)%in%c("Age_numeric")))==1 & length(table(na.omit(data$Age_numeric)))!=1){

fit1<- lm(trait ~ Age_numeric, data=data)
fit2<- lm(trait ~ Age_numeric+I(Age_numeric^2), data=data)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age_numeric"]<0.05){
data$trait<-resid(fit1)}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age_numeric^2)"]<0.05){
data$trait<-resid(fit2)}
}
#standardise
nmiss<-which(!is.na(data[,"trait"]))
data[nmiss,"trait"]<-(data[nmiss,"trait"]-mean(data[nmiss,"trait"]))/sd(data[nmiss,"trait"])

outdata<-data.frame(IID=data$IID,trait=data$trait)
m<-match(IID[,2],outdata$IID)
outdata<-data.frame(IID=IID[,2],trait=outdata$trait[m])

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

male <- data[which(data$Sex_factor=="M"),]
female <- data[which(data$Sex_factor=="F"),]

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

#adjust for covariates

if(length(which(names(data)%in%"Age_numeric"))>0 & length(table(na.omit(data$Age_numeric)))!=1){

fit1<- lm(trait ~ Age_numeric, data=male)
fit2<- lm(trait ~ Age_numeric+I(Age_numeric^2), data=male)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age_numeric"]<0.05){
male$trait<-resid(fit1)}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age_numeric^2)"]<0.05){
male$trait<-resid(fit2)}

fit1<- lm(trait ~ Age_numeric, data=female)
fit2<- lm(trait ~ Age_numeric+I(Age_numeric^2), data=female)

if(coefficients(summary(fit1))[,"Pr(>|t|)"]["Age_numeric"]<0.05){
female$trait<-resid(fit1)}

if(coefficients(summary(fit2))[,"Pr(>|t|)"]["I(Age_numeric^2)"]<0.05){
female$trait<-resid(fit2)}

}

#standardise
nmiss_male<-which(!is.na(male[,"trait"]))
male[nmiss_male,"trait"]<-(male[nmiss_male,"trait"]-mean(male[nmiss_male,"trait"]))/sd(male[nmiss_male,"trait"])
nmiss_female<-which(!is.na(female[,"trait"]))
female[nmiss_female,"trait"]<-(female[nmiss_female,"trait"]-mean(female[nmiss_female,"trait"]))/sd(female[nmiss_female,"trait"])

outdata<-rbind(male,female)
m<-match(IID[,2],outdata$IID)
outdata<-data.frame(IID=IID[,2],trait=outdata$trait[m])

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
}
null <- dev.off()

colnames(outdata.all)<-c("IID",traits)
m<-match(IID[,2],outdata.all[,1])
outdata.all<-data.frame(FID=IID[,1],outdata.all[m,])
write.table(outdata.all[,-1],file=transformed.phenotypes,sep="\t",quote=F,row.names=F,col.names=T)
write.table(outdata.all,file=phenotypes.plink,quote=F,row.names=F,col.names=F)
