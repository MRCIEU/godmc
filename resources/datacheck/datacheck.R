errorlist <- list()
warninglist <- list()


message("Checking R version")
currentr <- paste0(R.Version()['major'], ".", R.Version()['minor'])
ch <- compareVersion(currentr, "3.1")
if(ch == -1)
{
	stop("You are running R version ", currentr, ". Please upgrade to at least 3.1.0.")
}



message("Checking that all required packages are present")

pkglist <- c(
	"lattice",
	"ggplot2",
	"data.table",
	"MatrixEQTL",
	"parallel",
	"GenABEL",
	"matrixStats",
	"plyr",
	"SNPRelate",
	"GENESIS",
	"meffil"
)

index <- pkglist %in% rownames(installed.packages())
if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
} else {
	message("All required packages installed")
}

# Check R version


library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
bim_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
betas_file <- as.character(args[3]);
covariates_file <- as.character(args[4]);
cellcounts_file <- as.character(args[5]);
ids <- as.character(args[6]);
ids_plink <- as.character(args[7]);
snpsbychr_file <- as.character(args[8]);
snpsbychr_plot <- as.character(args[9]);
controlsnps_file <- as.character(args[10]);
phenotypes_file <- as.character(args[11]);
cnv_file <- as.character(args[12]);
cohort_descriptives_file <- as.character(args[13])
age_distribution_plot <- as.character(args[14])
ewas_phenotype_list_file <- as.character(args[15])
quality_file <- as.character(args[16])
quality_plot <- as.character(args[17])
methylation_summary_file <- as.character(args[18])

# BIM file check
message("Checking bim file: ", bim_file)
controlsnps <- read.table(controlsnps_file, header=F, stringsAsFactors=F)
bim <- as.data.frame(fread(bim_file))

message("Number of SNPs: ", nrow(bim))

# test chr coding
chrno <- table(bim[,1])
w <- which(! names(chrno) %in% as.character(c(1:22)))

print(data.frame(chrno))

if(length(w) > 0)
{
	msg <- paste0("There are some chromosomes other than 1-22, they will be removed")
	warninglist <- c(warninglist, msg)
	message("Warning: ", msg)
}

w <- which(names(chrno) %in% c(1:22))
if(length(w)<22)
{
	msg <- "Please change chromosome coding to 1-22, please dont use chr1, chr2 etc."
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


message("Checking strand")

bim2<-data.frame(bim,alleles=paste(bim[,5],bim[,6],sep=""))
w<-which(bim2$alleles=="GA")
bim2$alleles[w]<-"AG"
w<-which(bim2$alleles=="CA")
bim2$alleles[w]<-"AC"
w<-which(bim2$alleles=="TC")
bim2$alleles[w]<-"CT"
w<-which(bim2$alleles=="TG")
bim2$alleles[w]<-"GT"

message("Checking strand against control SNPs")
for (i in 1:22)
{
	chr <- bim2[which(bim2[,1] %in% i),]
	controlsnps.chr <- na.omit(controlsnps[which(controlsnps$V2 %in% i), ])
	controlsnps.chr$V6<-as.character(controlsnps.chr$V6)
	w<-which(controlsnps.chr$V5=="-"&controlsnps.chr$V6=="AG")
	if(length(w)>0){controlsnps.chr$V6[w]<-"CT"}
    w<-which(controlsnps.chr$V5=="-"&controlsnps.chr$V6=="AC")
	if(length(w)>0){controlsnps.chr$V6[w]<-"GT"}

	m<-match(controlsnps.chr$V3,chr$V4)
	chr<-chr[na.omit(m),]
	m<-match(chr$V4,controlsnps.chr$V3)
	controlsnps.chr<-controlsnps.chr[m,]
    strand.check<-sum(controlsnps.chr$V6==chr$alleles,na.rm=T)/nrow(controlsnps.chr)
    message("Chr ", i, " proportion in agreement: ", strand.check)	
	if(strand.check<0.95)
	{
		msg <- paste0("please check strand for chromosome ",i," as more than 5% of your SNPs have strand issues")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
}


message("Checking for duplicate SNPs")
if(any(duplicated(bim[,2])))
{
	msg <- "duplicate SNPs in bim file"
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

###
# id <- paste(bim[,1], bim[,4], nchar(as.character(bim[,5])),nchar(as.character(bim[,6])))
# l1 <- length(id)
# l2 <- length(unique(id))
# l <- (length(id)-length(unique(id)))/length(id)
# if(l > 0)
# {
# 	warning("WARNING: you have duplicated CHR:POS:{SNP/INDEL} positions. See the wiki for info on how to handle this.")
# }

#######################################################

message("Checking build against control SNPs by checking positions of common controlSNPs")
no.SNPs.bychr <- NULL
for (i in 1:22)
{
	chr <- bim[which(bim[,1] %in% i),]
	no.SNPs <- nrow(chr)
	controlsnps.chr <- controlsnps[which(controlsnps$V2 %in% i), ]
	w <- which(as.character(chr[,4]) %in% controlsnps.chr$V3)
	pos.check <- length(w)/nrow(controlsnps.chr)
	message("Chr ", i, " proportion in agreement: ", pos.check)
	no.SNPs.bychr <- append(no.SNPs.bychr, no.SNPs)
	if(pos.check<0.50)
	{
		msg <- paste0("please change positions for chromosome ",i, " to build 37 as less than 50% of common controlsnps are found")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
}
pdf(snpsbychr_plot, height=6, width=6)
barplot(no.SNPs.bychr, main="no of SNVs by chromosome",xlab="chromosome",names=c(1:22),cex.names=0.6,cex.axis=0.6)
null <- dev.off()

write.table(no.SNPs.bychr,snpsbychr_file,sep="\t",quote=F,row.names=F,col.names=F)



##

message("Checking imputation quality scores: ", quality_file)
qual <- as.data.frame(fread(quality_file,header=T))

if(ncol(qual) != 3)
{
	msg <- paste0("Expecting 3 columns in the imputation quality file: SNP ID, MAF and quality score.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(any(qual[,2] < 0) | any(qual[,2] > 1))
{
	msg <- paste0("Second column of quality scores file should be MAF. Some of the provided values fall outside the range of 0-1")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(any(qual[,3] > 1.1))
{
	msg <- paste0("third column of quality scores file should be the info score. Some of the provided values are above 1.")
	msg <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

prop <- sum(bim[,2] %in% qual[,1]) / nrow(bim)
if(prop < 0.95)
{
	msg <- paste0("Less then 95% of SNPs in the genetic data have info scores provided.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

message(round(prop*100, 2), "% of the SNPs in the data have matching info scores.")
qual <- qual[qual[,1] %in% bim[,2], ]
message("Retaining ", nrow(qual), " quality scores.")


names(qual) <- c("V1", "V2", "V3")
index <- qual$V2 > 0.5
qual$V2[index] <- 1 - qual$V2[index]

prop <- sum(qual[,2] < 0.01) / nrow(qual)
if(prop > 0.1)
{
	msg <- paste0("more than 10% of the retained quality scores have a MAF < 0.01. Please filter on MAF < 0.01")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

prop <- sum(qual[,3] < 0.8) / nrow(qual)
if(prop > 0.1)
{
	msg <- paste0("more than 10% of the retained quality scores have a quality score < 0.4. Please filter the data. The wiki has a guide for doing this.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

message("Plotting info scores")
qual$mafbin <- cut(qual[,2], breaks=30)
pdf(quality_plot)
boxplot(V3 ~ mafbin, qual, xlab="MAF", ylab="Imputation quality")
null <- dev.off()

##
#FAM file check

#Family ID
#Individual ID
#Paternal ID
#Maternal ID
#Sex (1=male; 2=female; other=unknown)
#Phenotype

message("Checking fam file: ", fam_file)

fam <- read.table(fam_file,header=F,stringsAsFactors=F)

if(any(duplicated(fam[,2])))
{
	msg <- paste0("Individual identifier is not unique. Please fix this before going on.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(any(grepl("_",fam[,2])))
{
	msg <- paste0("please remove underscores from individual ids")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


#METHYLATION data check
#indiv1 indiv2
#CpG1
#CpG2

message("Checking methylation data: ", betas_file)
load(betas_file)

if(! "norm.beta" %in% ls())
{
	msg <- paste0("please save methylation matrix with the object name norm.beta")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("norm.beta object found")
#is methylation data a matrix

if(!is.matrix(norm.beta))
{
	msg <- paste0("please transform methylation norm.beta to a matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
#are CpGs in rows?
d1 <- nrow(norm.beta)
nid_meth <- ncol(norm.beta)
message("Number of individuals with methylation data: ", nid_meth)
message("Number of CpGs: ", d1)
if(d1 < nid_meth)
{
	msg <- paste0("please transpose methylation matrix (CpGs in rows; samples in columns)")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("Data is a correctly oriented matrix")

#are individuals unique
if(any(duplicated(colnames(norm.beta))))
{
	msg <- paste0("please remove duplicate samples from methylation data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("No duplicate IDs")

#check for NAs in beta matrix
if(any(is.na(norm.beta)))
{
	msg <- paste0("please remove NAs from methylation matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("No NAs in data")

#check for negative values in beta matrix
if(any(norm.beta < 0))
{
	msg <- paste0("please remove negative values from methylation matrix. Are these beta values?")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#check for values above 1 in beta matrix
if(any(norm.beta > 1)) 
{
	msg <- paste0("please remove values > 1 from methylation matrix. Are these beta values?")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("All values are within 0-1")

if(any(grepl("rs", rownames(norm.beta))))
{
	msg <- paste0("there are SNPs in the methylation data. Please remove all rows with rs IDs.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#extract list of individuals with geno+methylation data
overlap <- intersect(colnames(norm.beta),fam[,2])
n.overlap <- length(overlap)
if(n.overlap < 50)
{
	msg <- paste0("fewer than 50 subjects with methylation and genotype data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

w <- which(fam[,2] %in% overlap)
fam2 <- fam[w,1:2]
message(nrow(fam2), " individuals present in both genetic and methylation datasets")
write.table(fam2[,2],ids,sep="\t",quote=F,row.names=F,col.names=F)
write.table(fam2[,1:2],ids_plink,sep="\t",quote=F,row.names=F,col.names=F)


#CELLCOUNTS

message("Checking cell counts data: ", cellcounts_file)
if(cellcounts_file != "NULL")
{
	cc <- read.table(cellcounts_file,header=T)
	c1 <- dim(cc)[1]
	c2 <- dim(cc)[2]

	if(c1!=nid_meth)
	{
		msg <- paste0("number of samples in cell counts file is not the same as in beta matrix")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)   
	}

	w <- which(names(cc)[1] %in% c("IID"))
	if(w!=1)
	{
		msg <- paste0("first column from cellcounts file should be the sample identifier with the name IID")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	if(c2<3)
	{
		msg <- paste0("are there any columns with cell counts missing in the cell counts file?")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	a <- apply(cc,2,function(x) y<-length(which(is.na(x))))
	if(length(which(a > 0.1*nid_meth)))
	{
		msg <- paste0("more than 10% of missingness in one of the cellcounts")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
	message("Number of cell types: ", c2)
	message("Cell types:\n", paste(names(cc)[-1], collapse="\n"))
} else {
	message("No cell counts are provided, these will be estimated by the pipeline.")
}

#COVARIATES
message("Checking covariates file: ", covariates_file)
covar <- read.table(covariates_file,header=T)
cov1 <- dim(covar)[1]
cov2 <- dim(covar)[2]
# if(cov1!=nid_meth)
# {
# 	warning("WARNING: number of samples in covariates file is not the same as in beta matrix")
# }

commonids_mgc <- Reduce(intersect, list(colnames(norm.beta), covar$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids_mgc))

if(length(commonids_mgc) < 50)
{
	msg <- paste0("must have at least 50 individuals with covariate, methylation and genetic data.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

w <- which(names(covar)[1] %in% c("IID"))
if(w!=1)
{
	msg <- paste0("first column from cellcounts file should be the sample identifier with the name IID")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(cov2<3)
{
	msg <- paste0("are there any covariates missing in the covariates file? Sex and Age are required")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

g1<-grep("_factor",names(covar))
g2<-grep("_numeric",names(covar))
g<-unique(c(g1,g2))

if(length(g)!=(cov2-1))
{
	msg <- paste0("have you specified whether your covariates are factors or numeric in the header of the covariates file? Please make sure your column headers are e.g. 'IID Sex_factor Age_numeric' etc")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


for (i in 1:length(g1))
{
	if(length(table(na.omit(covar[,g1[i]])))==(cov1))
	{
		msg <- paste0(g1[i], " is specified as a factor but has the same number of levels as individuals")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}
}

a <- apply(covar,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*nid_meth)))
{
	msg <- paste0("more than 10% of missingness in the following covariates:\n", 
		paste(names(covar)[a>0.1*nid_meth], collapse="\n")
	)
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(! "Sex_factor" %in% names(covar))
{
	msg <- paste0("There is no Sex_factor variable in the covariate file. Please provide M/F values, even if they are all the same sex.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(any(is.na(covar$Sex_factor)))
{
	msg <- paste0("There are some missing values in the Sex_factor column. Please make sure all individuals have data for this column.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

index <- covar$Sex_factor %in% c("M", "F")
if(any(!index))
{
	msg <- paste0("There are some values in the Sex_factor column that are neither M nor F. Please make sure all individuals have data for this column.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(! "Age_numeric" %in% names(covar))
{
	msg <- paste0("There is no Age_numeric variable in the covariate file. Please provide age in years, even if they are all the same age.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

pdf(age_distribution_plot, height=6, width=6)
hist(covar$Age_numeric, breaks=50, xlab="Age", main=paste("age distribution (N=", length(which(!is.na(covar$Age_numeric))),")",sep=""),cex.main=0.7)
null <- dev.off()


if(any(is.na(covar$Age_numeric)))
{
	msg <- paste0("Some individuals don't have ages. Please make sure there are no missing values.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(any(covar$Age_numeric < 0))
{
	msg <- paste0("Some negative values in the age column.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(mean(covar$Age_numeric, na.rm=T) > 100)
{
	msg <- paste0("Average age is above 100, please make sure age is provided in years.")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("Number of individuals with covariate data: ", cov1)
message("Covariates provided:\n", paste(names(covar)[-1], collapse="\n"))
message("Average age: ", mean(covar$Age_numeric))
message("Number of males: ", sum(covar$Sex_factor == "M"))
message("Number of females: ", sum(covar$Sex_factor == "F"))


#EWAS phenotypes

message("Checking phenotypes: ", phenotypes_file)
if(phenotypes_file != "NULL")
{
	ph<-read.table(phenotypes_file,header=T)
	p1<-dim(ph)[1]
	p2<-dim(ph)[2]

	# if(p1!=nid_meth)
	# {
	# 	warning("WARNING: number of samples in phenotype file is not the same as in beta matrix","\n")   
	# }
	
	if(names(ph)[1] != "IID")
	{
		msg <- paste0("first column from phenotype file should be the sample identifier with the name IID")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	commonids_mpc <- Reduce(intersect, list(covar$IID, colnames(norm.beta), ph$IID))
	message(length(commonids_mpc), " in common between covariate, methylation and phenotype data")

	if(length(commonids_mpc) < 50)
	{
		msg <- paste0("fewer than 50 subjects with methylation, CNV and covariate data")
		warninglist <- c(warninglist, msg)
		warning("Warning: ", msg)
	}

	ph <- subset(ph, IID %in% commonids_mpc)

	if(p2 < 2)
	{
		msg <- paste0("No phenotypes present. Please set the phenotype variable to 'NULL' in the config file")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}

	nom <- names(ph)[-1][names(ph)[-1] %in% c("BMI", "Height")]
	if(length(nom) < 1)
	{
		msg <- paste0("Neither 'Height' nor 'BMI' variables are present in the phenotype file. Please check that the columns are correctly entered (note capitalisation).")
		errorlist <- c(errorlist, msg)
		warning("ERROR: ", msg)
	}


	a <- apply(ph,2,function(x) sum(is.na(x), na.rm=T))
	if(any(a > 0.1*nid_meth))
	{
		nom <- names(ph)[a > 0.1*nid_meth]
		msg <- paste0("more than 10% of missingness in the following EWAS phenotypes:\n", paste(nom, collapse="\n"))
		warninglist <- c(warninglist, msg)
		warning("Warning: ", msg)
	}

	if("Height" %in% nom)
	{
		message("Checking Height")
		m1 <- mean(ph$Height,na.rm=T)
		if(m1<1.0|m1>2.5)
		{
			msg <- paste0("please convert Height units to metres")
			errorlist <- c(errorlist, msg)
			warning("ERROR: ", msg)
		}
	}

	if("BMI" %in% nom)
	{
		message("Checking BMI")
		m1<-mean(ph$BMI,na.rm=T)
		if(m1<10|m1>35)
		{
			msg <- paste0("please convert BMI units to kg/m2")
			errorlist <- c(errorlist, msg)
			warning("ERROR: ", msg)
		}
	}

	write.table(names(ph)[-1], file=ewas_phenotype_list_file, row=F, col=F, qu=F)
} else {
	msg <- paste0("No phenotypes have been provided.\nWARNING: EWAS will not be performed.")
	warninglist <- c(warninglist, msg)
	message("WARNING: ", msg)
}

#CNV data check
#indiv1 indiv2
#cnv1
#cnv2


message("Checking CNV data: ", cnv_file)
load(cnv_file)

if(! "cnv" %in% ls())
{
	msg <- paste0("please save cnv matrix with the object name 'cnv'")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#is cnv data a matrix

if(!is.matrix(cnv))
{
	msg <- paste0("please transform cnv object to a matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#are CpGs in rows?
d1 <- dim(cnv)[1]
d2 <- dim(cnv)[2]
if(d1 < d2)
{
	msg <- paste0("please transpose cnv matrix (cnvs in rows; samples in columns)")
	errorlist <- c(errorlist, msg)
    warning("ERROR: ", msg)
}

message("Number of positions with CNV data: ", d1)
message("Number of individuals in CNV data: ", d2)

#are individuals unique
if(any(duplicated(colnames(cnv))))
{
	msg <- paste0("please remove duplicate samples from cnv data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

#check for NAs in beta matrix
if(any(is.na(cnv)))
{
	msg <- paste0("please remove NAs from cnv matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

commonids_mcc <- Reduce(intersect, list(colnames(cnv), covar$IID, colnames(norm.beta)))
message(length(commonids_mcc), " samples in common between CNV, covariate and methylation data")

if(length(commonids_mcc) < 50)
{
	msg <- paste0("fewer than 50 subjects with methylation, CNV and covariate data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(d2 != dim(norm.beta)[2])
{
	msg <- paste0("number of subject with CNV data is not the same as number of subjects with methylation data")
	errorlist <- c(errorlist, msg)
    warning("ERROR: ", msg)
}

# Cohort characteristics

# Only count those that will be used in the analysis
covar <- subset(covar, IID %in% commonids_mgc)

cohort_summary <- list()
cohort_summary$mqtl_sample_size <- length(commonids_mgc)
cohort_summary$mcnv_sample_size <- length(commonids_mcc)
cohort_summary$n_males <- sum(covar$Sex_factor == "M",na.rm=T)
cohort_summary$n_females <- sum(covar$Sex_factor == "F",na.rm=T)
cohort_summary$mean_age <- mean(covar$Age_numeric,na.rm=T)
cohort_summary$median_age <- median(covar$Age_numeric,na.rm=T)
cohort_summary$sd_age <- sd(covar$Age_numeric,na.rm=T)
cohort_summary$max_age <- max(covar$Age_numeric,na.rm=T)
cohort_summary$min_age <- min(covar$Age_numeric,na.rm=T)
cohort_summary$Height_sample_size <- sum(!is.na(ph$Height))
cohort_summary$mean_Height <- mean(ph$Height,na.rm=T)
cohort_summary$median_Height <- median(ph$Height,na.rm=T)
cohort_summary$sd_Height <- sd(ph$Height,na.rm=T)
cohort_summary$max_Height <- max(ph$Height,na.rm=T)
cohort_summary$min_Height <- min(ph$Height,na.rm=T)
cohort_summary$BMI_sample_size <- sum(!is.na(ph$BMI))
cohort_summary$mean_BMI <- mean(ph$BMI,na.rm=T)
cohort_summary$median_BMI <- median(ph$BMI,na.rm=T)
cohort_summary$sd_BMI <- sd(ph$BMI,na.rm=T)
cohort_summary$max_BMI <- max(ph$BMI,na.rm=T)
cohort_summary$min_BMI <- min(ph$BMI,na.rm=T)
cohort_summary$n_CpGs <- nrow(norm.beta)
cohort_summary$n_SNP <- nrow(bim)
cohort_summary$covariates <- names(covar)[-1]


summariseMeth <- function(X, outlier_threshold, niter)
{
	message("Counting outliers in methylation matrix")
	
    norm.beta.copy <- X
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta.copy, na.rm=T)
		means <- rowMeans(norm.beta.copy, na.rm=T)
		norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	}
	outliers <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	
	message("Estimating means")
	means <- rowMeans(norm.beta.copy, na.rm=T)

	message("Estimating SDs")
	sds <- rowSds(norm.beta.copy, na.rm=T)

	message("Estimating medians")
	medians <- rowMedians(norm.beta.copy, na.rm=T)

	dat <- data.frame(cpg=rownames(norm.beta.copy), mean=means, median=medians, sd=sds, outlier=outliers)
	return(dat)
}

message("Generating summary stats of methylation")

meth_summary <- summariseMeth(norm.beta, 10, 3)

save(cohort_summary, file=cohort_descriptives_file)
save(meth_summary, file=methylation_summary_file)

message("\n\nCompleted checks\n")

message("Summary of data:")
for(i in 1:length(cohort_summary))
{
	a <- cohort_summary[[i]]
	if(is.numeric(a)) a <- round(a, 2)
	message(names(cohort_summary)[i], ": ", paste(a, collapse=", "))
}


if(length(warninglist) > 0)
{
	message("\n\nPlease take note of the following warnings, and fix and re-run the data check if you see fit:")
	null <- sapply(warninglist, function(x)
	{
		message("- ", x)
	})
}


if(length(errorlist) > 0)
{
	message("\n\nThe following errors were encountered, and must be addressed before continuing:")
	null <- sapply(errorlist, function(x)
	{
		message("- ", x)
	})
	q(status=1)
}
message("\n\n")