args = (commandArgs(TRUE));
bim= as.character(args[1]);
fam= as.character(args[2]);
betas= as.character(args[3]);
covariates=as.character(args[4]);
cellcounts=as.character(args[5]);
ids=as.character(args[6]);
ids_plink=as.character(args[7]);
SNPsbychr=as.character(args[8]);
no.SNPs.by.chr.plot=as.character(args[9]);
snpsforpositioncheck=as.character(args[10]);
EWASphenotypes=as.character(args[11]);
height=as.character(args[12]);
bmi=as.character(args[13]);
cnvs=as.character(args[14]);
log_directory=as.character(args[15])

#BIM file check
controlsnps<-read.table(snpsforpositioncheck,header=F,stringsAsFactors=F)
bim<-read.table(bim,header=F,stringsAsFactors=F)

#CHR
#ID
#cM
#BP
#A1
#A2

#test chr coding
chrno<-table(bim[,1])
w<-which(names(chrno)%in%as.character(c(1:22))==F)
cat("check chromosome coding","\n")

if(length(w)>0){
message("Warning: There are some chromosomes other than 1-22, they will be removed","\n")
}

w<-which(names(chrno)%in%c(1:22)==T)
if(length(w)<22) {
stop("ERROR:please change chromosome coding to 1-22","\n")
}

a1<-data.frame(table(bim[,5]))
a2<-data.frame(table(bim[,6]))

allele.out<-NULL
for (t in 1:dim(a1)[1]){
testa1<-unlist(strsplit(as.character(a1$Var1[t]),split="",fixed=T))
w<-which(testa1%in%c("A","C","T","G")==F)
if(length(w)>0){
allele<-which(bim[,5]%in%a1$Var[t])
allele.out<-rbind(allele.out,bim[allele,])}}

if(dim(allele.out)[1]>0){perc.miscoding<-dim(allele.out)[1]/dim(bim)[1]}
if(perc.miscoding>0.01){stop("ERROR:more than 1% of miscoding alleles, please change allele coding to A,C,T,G","\n")}


allele.out<-NULL
for (t in 1:dim(a1)[1]){
testa2<-unlist(strsplit(as.character(a2$Var1[t]),split="",fixed=T))
w<-which(testa2%in%c("A","C","T","G")==F)
if(length(w)>0){
allele<-which(bim[,6]%in%a2$Var[t])
allele.out<-rbind(allele.out,bim[allele,])}}

if(dim(allele.out)[1]>0){perc.miscoding<-dim(allele.out)[1]/dim(bim)[1]}
if(perc.miscoding>0.01){stop("ERROR:more than 1% of miscoding alleles, please change allele coding to A,C,T,G","\n")}

###
id<-paste(bim[,1],bim[,4],nchar(as.character(bim[,5])),nchar(as.character(bim[,6])))
l1<-length(id)
l2<-length(unique(id))
l<-(length(id)-length(unique(id)))/length(id)
if(l>0){
stop("ERROR:you have duplicated CHR:POS:{SNP/INDEL} positions. See the wiki for info on how to handle this.","\n")  
}

#######################################################

no.SNPs.bychr<-NULL
for (i in 1:22){
cat("check position and alleles for chromosome ",i,"\n") 
chr<-bim[which(bim[,1]%in%i),]
no.SNPs<-dim(chr)[1]
controlsnps.chr<-controlsnps[which(controlsnps$V2%in%i),]
w<-which(chr[,4]%in%controlsnps.chr$V3)
pos.check<-length(w)/dim(controlsnps.chr)[1]
cat(pos.check,"\n")
no.SNPs.bychr<-append(no.SNPs.bychr,no.SNPs)
if(pos.check<0.80){
stop("ERROR:please change positions for chromosome ",i,"\n")   
}
}
pdf(paste(no.SNPs.by.chr.plot,".pdf",sep=""),height=6,width=6)
barplot(no.SNPs.bychr, main="no of SNVs by chromosome",xlab="chromosome",names=c(1:22),cex.names=0.6,cex.axis=0.6)
dev.off()

write.table(no.SNPs.bychr,SNPsbychr,sep="\t",quote=F,row.names=F,col.names=F)


##
#FAM file check

#Family ID
#Individual ID
#Paternal ID
#Maternal ID
#Sex (1=male; 2=female; other=unknown)
#Phenotype

fam<-read.table(fam,header=F,stringsAsFactors=F)

d1<-length(fam[,2])
d2<-length(unique(fam[,2]))
if(d1!=d2){
stop("ERROR:individual identifier is not unique","\n")   
}

g<-grep("_",fam[,2])
if(length(g)>0){
stop("ERROR:please remove underscores from individual ids","\n")   
}


#METHYLATION data check
#indiv1 indiv2
#CpG1
#CpG2



load(betas)

l<-ls()
w<-which(l%in%c("norm.beta"))

if(length(w)<1){
stop("ERROR:please save methylation matrix with the name norm.beta","\n")
}

#is methylation data a matrix

if(is.matrix(norm.beta)==F){
stop("ERROR:please transform methylation norm.beta to a matrix","\n")	
}

#are CpGs in rows?
d1<-dim(norm.beta)[1]
d2<-dim(norm.beta)[2]
if(d1<d2){
stop("ERROR:please transpose methylation matrix (CpGs in rows; samples in columns)","\n")   
}

#are individuals unique
c1<-length(colnames(norm.beta))
c2<-length(unique(colnames(norm.beta)))

if(c1>c2){stop("ERROR:please remove duplicate samples from methylation data")}

#check for NAs in beta matrix
if(any(is.na(norm.beta))){stop("ERROR:please remove NAs from methylation matrix","\n")}

#check for negative values in beta matrix
if(any(norm.beta < 0)) {stop("ERROR:please remove negative values from methylation matrix","\n")}

#check for values above 1 in beta matrix
if(any(norm.beta > 1)) {stop("ERROR:please remove negative values from methylation matrix","\n")}

#extract list of individuals with geno+methylation data
overlap<-intersect(colnames(norm.beta),fam[,2])
n.overlap<-length(overlap)
if(n.overlap<100){
stop("ERROR:less than 100 subjects with methylation and genotype data","\n")
}
w<-which(fam[,2]%in%overlap)
fam2<-fam[w,1:2]
write.table(fam2[,2],ids,sep="\t",quote=F,row.names=F,col.names=F)
write.table(fam2[,1:2],ids_plink,sep="\t",quote=F,row.names=F,col.names=F)


#CELLCOUNTS
if(cellcounts != "NULL")){
cc<-read.table(cellcounts,header=T)
c1<-dim(cc)[1]
c2<-dim(cc)[2]
if(c1!=d2){
stop("ERROR:number of samples in cell counts file is not the same as in beta matrix","\n")   
} else {
	message("No cell counts are provided, these will be estimated by the pipeline.")
}

w<-which(names(cc)[1]%in%c("IID"))
if(w!=1){
stop("ERROR:first column from cellcounts file should be the sample identifier with the name IID","\n")
}

if(c2<3){
stop("ERROR:are there any columns with cell counts missing in the cell counts file?","\n")   
}

a<-apply(cc,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*d2))){
stop("ERROR:more than 10% of missingness in one of the cellcounts","\n")
}
}

#COVARIATES
covar<-read.table(covariates,header=T)
cov1<-dim(covar)[1]
cov2<-dim(covar)[2]
if(cov1!=d2){
stop("ERROR:number of samples in covariates file is not the same as in beta matrix","\n")   
}

w<-which(names(covar)[1]%in%c("IID"))
if(w!=1){
stop("ERROR:first column from cellcounts file should be the sample identifier with the name IID","\n")
}

if(cov2<3){
stop("ERROR:are there any covariates missing in the covariates file?","\n")   
}

a<-apply(covar,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*d2))){
stop("ERROR:more than 10% of missingness in one of the covariates","\n")   
}

if(! "Sex" %in% names(covar))
{
	stop("ERROR: There is no Sex variable in the covariate file. Please provide M/F values, even if they are all the same sex.")
}

if(any(is.na(covar$Sex)))
{
	stop("ERROR: There are some values in the Sex column that are neither M nor F. Please make sure all individuals have data for this column.")
}

index <- covar$Sex == "M" | covar$sex == "F"
if(any(!index))
{
	stop("ERROR: There are some values in the Sex column that are neither M nor F. Please make sure all individuals have data for this column.")
}

if(! "Age" %in% names(covar))
{
	stop("ERROR: There is no Age variable in the covariate file. Please provide age in years, even if they are all the same age.")
}

if(any(is.na(covar$Age)))
{
	stop("ERROR: Some individuals don't have ages. Please make sure there are no missing values.")
}

if(any(covar$Age < 0))
{
	stop("ERROR: Some negative values in the age column.")
}

if(mean(covar$Age, na.rm=T) > 100)
{
	stop("ERROR: Average age is above 100, please make sure age is provided in years.")
}




#EWAS phenotypes
if(!is.na(EWASphenotypes)){
ph<-read.table(EWASphenotypes,header=T)
p1<-dim(ph)[1]
p2<-dim(ph)[2]
if(p1!=d2){
stop("ERROR:number of samples in phenotype file is not the same as in beta matrix","\n")   
}

w<-which(names(ph)[1]%in%c("IID"))
if(w!=1){
stop("ERROR:first column from phenotype file should be the sample identifier with the name IID","\n")
}

if(p2<2){
stop("ERROR:are there any columns with phenotypes missing in the EWAS phenotypes file?","\n")   
}

a<-apply(ph,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*d2))){
stop("ERROR:more than 10% of missingness in one of the EWASphenotypes","\n")
}

if(height != "no")){
w<-which(names(ph)%in%c("Height"))
if(length(w)<1){
    stop("ERROR:please change height variable to be called 'Height' (note capitalisation)","\n")
}
if(length(w)>0){
m1<-mean(ph[,w],na.rm=T)
if(m1<1.0|m1>2.5){
stop("ERROR:please convert Height units to centimeters","\n")
}
}
}

if(bmi != "no")){
w<-which(names(ph)%in%c("BMI"))
if(length(w)<1){
    stop("ERROR:please change BMI variable to 'BMI' (note capitalisation)","\n")
}

if(length(w)>0){
m1<-mean(ph[,w],na.rm=T)
if(m1<10|m1>35){
stop("ERROR:please convert BMI units to kg/m2","\n")
}
}
}

a<-apply(ph,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*d2))){
stop("ERROR:more than 10% of missingness in one of the EWASphenotypes","\n")
}
}

#CNV data check
#indiv1 indiv2
#cnv1
#cnv2



load(cnvs)

l<-ls()
w<-which(l%in%c("cnv"))

if(length(w)<1){
    stop("ERROR:please save cnv matrix with the name cnv","\n")
}

#is cnv data a matrix

if(is.matrix(cnv)==F){
    stop("ERROR:please transform cnv object to a matrix","\n")
}

#are CpGs in rows?
d1<-dim(cnv)[1]
d2<-dim(cnv)[2]
if(d1<d2){
    stop("ERROR:please transpose cnv matrix (cnvs in rows; samples in columns)","\n")
}

#are individuals unique
c1<-length(colnames(cnv))
c2<-length(unique(colnames(cnv)))

if(c1>c2){stop("ERROR:please remove duplicate samples from cnv data")}

#check for NAs in beta matrix
if(any(is.na(cnv))){stop("ERROR:please remove NAs from cnv matrix","\n")}


# Cohort characteristics

cohort_summary <- list()
cohort_summary$sample_size <- length(ids)
cohort_summary$n_males <- sum(covar$Sex == "M")
cohort_summary$n_females <- sum(covar$Sex == "F")
cohort_summary$mean_age <- mean(covar$Age)
cohort_summary$sd_age <- sd(covar$Age)
cohort_summary$max_age <- max(covar$Age)
cohort_summary$min_age <- min(covar$Age)
cohort_summary$mean_Height <- mean(ph$Height)
cohort_summary$sd_Height <- sd(ph$Height)
cohort_summary$max_Height <- max(ph$dHeight)
cohort_summary$min_Height <- min(ph$Height)
cohort_summary$mean_BMI <- mean(ph$BMI)
cohort_summary$sd_BMI <- sd(ph$BMI)
cohort_summary$max_BMI <- max(ph$dBMI)
cohort_summary$min_BMI <- min(ph$BMI)
cohort_summary$n_snp <- nrow(bim)

save(cohort_summary, file=file.path(log_directory, "cohort_descriptives.RData"))

cat("You successfully performed all datachecks!","\n")







