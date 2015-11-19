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
	stop("The following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
} else {
	message("All required packages installed")
}


library(data.table)

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
# bim <- read.table(bim_file, stringsAsFactors=F, colClass=c("numeric", "character", "numeric", "numeric", "character", "character"))
bim <- as.data.frame(fread(bim_file))

message("Number of SNPs: ", nrow(bim))

# test chr coding
chrno <- table(bim[,1])
w <- which(! names(chrno) %in% as.character(c(1:22)))

print(data.frame(chrno))

if(length(w) > 0)
{
	message("Warning: There are some chromosomes other than 1-22, they will be removed","\n")
}

w <- which(names(chrno) %in% c(1:22))
if(length(w)<22)
{
	stop("ERROR: please change chromosome coding to 1-22")
}

message("Checking alleles")
a1 <- data.frame(table(bim[,5]))
a2 <- data.frame(table(bim[,6]))

allele.out <- NULL
for (t in 1:nrow(a1))
{
	testa1 <- unlist(strsplit(as.character(a1$Var1[t]), split="", fixed=T))
	w <- which(! testa1 %in% c("A","C","T","G"))
	if(length(w) > 0)
	{
		allele <- which(bim[,5] %in% a1$Var[t])
		allele.out <- rbind(allele.out,bim[allele,])
	}
}

if(!is.null(allele.out))
{
	perc.miscoding <- nrow(allele.out)/nrow(bim)
	if(perc.miscoding > 0.01)
	{
		stop("ERROR: more than 1% of miscoding alleles, please change allele coding to A,C,T,G")
	}
}

allele.out <- NULL
for (t in 1:nrow(a2))
{
	testa2 <- unlist(strsplit(as.character(a2$Var1[t]), split="", fixed=T))
	w <- which(! testa2 %in% c("A","C","T","G"))
	if(length(w) > 0)
	{
		allele <- which(bim[,5] %in% a2$Var[t])
		allele.out <- rbind(allele.out,bim[allele,])
	}
}

if(!is.null(allele.out))
{
	perc.miscoding <- nrow(allele.out)/nrow(bim)
	if(perc.miscoding > 0.01)
	{
		stop("ERROR: more than 1% of miscoding alleles, please change allele coding to A,C,T,G")
	}
}

message("Checking for duplicate SNPs")
if(any(duplicated(bim[,2])))
{
	stop("ERROR: duplicate SNPs in bim file")
}

###
# id <- paste(bim[,1], bim[,4], nchar(as.character(bim[,5])),nchar(as.character(bim[,6])))
# l1 <- length(id)
# l2 <- length(unique(id))
# l <- (length(id)-length(unique(id)))/length(id)
# if(l > 0)
# {
# 	stop("ERROR: you have duplicated CHR:POS:{SNP/INDEL} positions. See the wiki for info on how to handle this.")
# }

#######################################################

message("Checking position and alleles for chromosome against control SNPs")
no.SNPs.bychr <- NULL
for (i in 1:22)
{
	chr <- bim[which(bim[,1] %in% i),]
	no.SNPs <- nrow(chr)
	controlsnps.chr <- controlsnps[which(controlsnps$V2 %in% i), ]
	w <- which(chr[,4] %in% controlsnps.chr$V3)
	pos.check <- length(w)/nrow(controlsnps.chr)
	message("Chr ", i, " proportion in agreement: ", pos.check)
	no.SNPs.bychr <- append(no.SNPs.bychr, no.SNPs)
	if(pos.check<0.80)
	{
		stop("ERROR: please change positions for chromosome ",i, " to build 37")
	}
}
pdf(snpsbychr_plot, height=6, width=6)
barplot(no.SNPs.bychr, main="no of SNVs by chromosome",xlab="chromosome",names=c(1:22),cex.names=0.6,cex.axis=0.6)
tmp <- dev.off()

write.table(no.SNPs.bychr,snpsbychr_file,sep="\t",quote=F,row.names=F,col.names=F)



##

message("Checking imputation quality scores: ", quality_file)
qual <- as.data.frame(fread(quality_file,header=T))

if(ncol(qual) != 3)
{
	stop("ERROR: Expecting 3 columns in the imputation quality file: SNP ID, MAF and quality score.")
}

if(any(qual[,2] < 0) | any(qual[,2] > 1))
{
	stop("ERROR: second column of quality scores file should me MAF. Some of the provided values fall outside the range of 0-1")
}

if(any(qual[,3] > 1))
{
	stop("ERROR: third column of quality scores file should be the info score. Some of the provided values are above 1.")
}

prop <- sum(bim[,2] %in% qual[,1]) / nrow(bim)
if(prop < 0.95)
{
	stop("ERROR: Less then 95% of SNPs in the genetic data have info scores provided.")
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
	stop("ERROR: more than 10% of the retained quality scores have a MAF < 0.01. Please filter on MAF < 0.01")
}

prop <- sum(qual[,3] < 0.4) / nrow(qual)
if(prop > 0.1)
{
	stop("ERROR: more than 10% of the retained quality scores have a quality score < 0.4. Please filter the data. The wiki has a guide for doing this.")
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
	stop("ERROR: individual identifier is not unique")
}

if(any(grepl("_",fam[,2])))
{
	stop("ERROR: please remove underscores from individual ids")
}


#METHYLATION data check
#indiv1 indiv2
#CpG1
#CpG2

message("Checking methylation data: ", betas_file)
load(betas_file)

if(! "norm.beta" %in% ls())
{
	stop("ERROR: please save methylation matrix with the object name norm.beta")
}
message("norm.beta object found")
#is methylation data a matrix

if(!is.matrix(norm.beta))
{
	stop("ERROR: please transform methylation norm.beta to a matrix")
}
#are CpGs in rows?
d1 <- nrow(norm.beta)
nid_meth <- ncol(norm.beta)
message("Number of individuals with methylation data: ", nid_meth)
message("Number of CpGs: ", d1)
if(d1 < nid_meth)
{
	stop("ERROR: please transpose methylation matrix (CpGs in rows; samples in columns)")
}
message("Data is a correctly oriented matrix")

#are individuals unique
if(any(duplicated(colnames(norm.beta))))
{
	stop("ERROR: please remove duplicate samples from methylation data")
}
message("No duplicate IDs")

#check for NAs in beta matrix
if(any(is.na(norm.beta)))
{
	stop("ERROR: please remove NAs from methylation matrix")
}
message("No NAs in data")

#check for negative values in beta matrix
if(any(norm.beta < 0))
{
	stop("ERROR: please remove negative values from methylation matrix. Are these beta values?")
}

#check for values above 1 in beta matrix
if(any(norm.beta > 1)) 
{
	stop("ERROR: please remove values > 1 from methylation matrix. Are these beta values?")
}
message("All values are within 0-1")

if(any(grepl("rs", rownames(norm.beta))))
{
	stop("ERROR: there are SNPs in the methylation data. Please remove all rows with rs IDs.")
}

#extract list of individuals with geno+methylation data
overlap <- intersect(colnames(norm.beta),fam[,2])
n.overlap <- length(overlap)
if(n.overlap < 50)
{
	stop("ERROR: fewer than 50 subjects with methylation and genotype data")
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
		stop("ERROR: number of samples in cell counts file is not the same as in beta matrix")   
	}

	w <- which(names(cc)[1] %in% c("IID"))
	if(w!=1)
	{
		stop("ERROR: first column from cellcounts file should be the sample identifier with the name IID")
	}

	if(c2<3)
	{
		stop("ERROR: are there any columns with cell counts missing in the cell counts file?")
	}

	a <- apply(cc,2,function(x) y<-length(which(is.na(x))))
	if(length(which(a > 0.1*nid_meth)))
	{
		stop("ERROR: more than 10% of missingness in one of the cellcounts")
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
# 	stop("ERROR: number of samples in covariates file is not the same as in beta matrix")
# }

commonids <- Reduce(intersect, list(colnames(norm.beta), covar$IID, fam[,2]))
message("Number of samples with covariate, methylation and genetic data: ", length(commonids))

if(length(commonids) < 50)
{
	stop("ERROR: must have at least 50 individuals with covariate, methylation and genetic data.")
}

w <- which(names(covar)[1] %in% c("IID"))
if(w!=1)
{
	stop("ERROR: first column from cellcounts file should be the sample identifier with the name IID")
}

if(cov2<3)
{
	stop("ERROR: are there any covariates missing in the covariates file? Sex and Age are required")
}

a <- apply(covar,2,function(x) y<-length(which(is.na(x))))
if(length(which(a>0.1*nid_meth)))
{
	stop("ERROR: more than 10% of missingness in one of the covariates","\n")   
}

if(! "Sex" %in% names(covar))
{
	stop("ERROR:  There is no Sex variable in the covariate file. Please provide M/F values, even if they are all the same sex.")
}

if(any(is.na(covar$Sex)))
{
	stop("ERROR:  There are some values in the Sex column that are neither M nor F. Please make sure all individuals have data for this column.")
}

index <- covar$Sex %in% c("M", "F")
if(any(!index))
{
	stop("ERROR:  There are some values in the Sex column that are neither M nor F. Please make sure all individuals have data for this column.")
}

if(! "Age" %in% names(covar))
{
	stop("ERROR:  There is no Age variable in the covariate file. Please provide age in years, even if they are all the same age.")
}

pdf(age_distribution_plot, height=6, width=6)
hist(covar$Age, breaks=50, xlab="Age", main=paste("age distribution (N=", length(which(!is.na(covar$Age))),")",sep=""),cex.main=0.7)
dev.off()


if(any(is.na(covar$Age)))
{
	stop("ERROR:  Some individuals don't have ages. Please make sure there are no missing values.")
}

if(any(covar$Age < 0))
{
	stop("ERROR:  Some negative values in the age column.")
}

if(mean(covar$Age, na.rm=T) > 100)
{
	stop("ERROR:  Average age is above 100, please make sure age is provided in years.")
}
message("Number of individuals with covariate data: ", cov1)
message("Covariates provided:\n", paste(names(covar)[-1], collapse="\n"))
message("Average age: ", mean(covar$Age))
message("Number of males: ", sum(covar$Sex == "M"))
message("Number of females: ", sum(covar$Sex == "F"))


#EWAS phenotypes

message("Checking phenotypes: ", phenotypes_file)
if(phenotypes_file != "NULL")
{
	ph<-read.table(phenotypes_file,header=T)
	p1<-dim(ph)[1]
	p2<-dim(ph)[2]

	# if(p1!=nid_meth)
	# {
	# 	stop("ERROR: number of samples in phenotype file is not the same as in beta matrix","\n")   
	# }
	
	if(names(ph)[1] != "IID")
	{
		stop("ERROR: first column from phenotype file should be the sample identifier with the name IID")
	}

	commonids <- Reduce(intersect, list(covar$IID, colnames(norm.beta), ph$IID))
	message(length(commonids), " in common between covariate, methylation and phenotype data")

	if(p2 < 2)
	{
		stop("ERROR: No phenotypes present. Please set the phenotype variable to 'NULL' in the config file")
	}

	nom <- names(ph)[-1][names(ph)[-1] %in% c("BMI", "Height")]
	if(length(nom) < 1)
	{
		stop("ERROR: Neither 'Height' nor 'BMI' variables are present in the phenotype file. Please check that the columns are correctly entered (note capitalisation).")
	}


	a <- apply(ph,2,function(x) sum(is.na(x), na.rm=T))
	if(any(a > 0.1*nid_meth))
	{
		nom <- names(ph)[a > 0.1*nid_meth]
		stop("ERROR: more than 10% of missingness in at least one of the EWAS phenotypes:\n", paste(nom, collapse="\n"))
	}

	if("Height" %in% nom)
	{
		message("Checking Height")
		m1 <- mean(ph$Height,na.rm=T)
		if(m1<1.0|m1>2.5)
		{
			stop("ERROR: please convert Height units to metres")
		}
	}

	if("BMI" %in% nom)
	{
		message("Checking BMI")
		m1<-mean(ph$BMI,na.rm=T)
		if(m1<10|m1>35)
		{
			stop("ERROR: please convert BMI units to kg/m2")
		}
	}

	write.table(names(ph)[-1], file=ewas_phenotype_list_file, row=F, col=F, qu=F)
} else {
	message("WARNING: No phenotypes have been provided.\nWARNING: EWAS will not be performed.")
}

#CNV data check
#indiv1 indiv2
#cnv1
#cnv2


message("Checking CNV data: ", cnv_file)
load(cnv_file)

if(! "cnv" %in% ls())
{
	stop("ERROR: please save cnv matrix with the object name 'cnv'")
}

#is cnv data a matrix

if(!is.matrix(cnv))
{
	stop("ERROR: please transform cnv object to a matrix")
}

#are CpGs in rows?
d1 <- dim(cnv)[1]
d2 <- dim(cnv)[2]
if(d1 < d2)
{
    stop("ERROR: please transpose cnv matrix (cnvs in rows; samples in columns)")
}

message("Number of positions with CNV data: ", d1)
message("Number of individuals in CNV data: ", d2)

#are individuals unique
if(any(duplicated(colnames(cnv))))
{
	stop("ERROR: please remove duplicate samples from cnv data")
}

#check for NAs in beta matrix
if(any(is.na(cnv)))
{
	stop("ERROR: please remove NAs from cnv matrix","\n")
}

commonids <- Reduce(intersect, list(colnames(cnv), covar$IID, colnames(norm.beta)))
message(length(commonids), " samples in common between CNV, covariate and methylation data")

# Cohort characteristics

cohort_summary <- list()
cohort_summary$sample_size <- length(ids)
cohort_summary$n_males <- sum(covar$Sex == "M")
cohort_summary$n_females <- sum(covar$Sex == "F")
cohort_summary$mean_age <- mean(covar$Age)
cohort_summary$median_age <- median(covar$Age)
cohort_summary$sd_age <- sd(covar$Age)
cohort_summary$max_age <- max(covar$Age)
cohort_summary$min_age <- min(covar$Age)
cohort_summary$mean_Height <- mean(ph$Height)
cohort_summary$median_Height <- median(ph$Height)
cohort_summary$sd_Height <- sd(ph$Height)
cohort_summary$max_Height <- max(ph$Height)
cohort_summary$min_Height <- min(ph$Height)
cohort_summary$mean_BMI <- mean(ph$BMI)
cohort_summary$median_BMI <- median(ph$BMI)
cohort_summary$sd_BMI <- sd(ph$BMI)
cohort_summary$max_BMI <- max(ph$BMI)
cohort_summary$min_BMI <- min(ph$BMI)
cohort_summary$n_snp <- nrow(bim)
cohort_summary$covariates <- names(covar)[-1]


summariseMeth <- function(X, outlier_threshold)
{
	require(matrixStats)

	message("Removing outliers")
	sds <- rowSds(X, na.rm=T)
	means <- rowMeans(X, na.rm=T)
	X[X > means + sds*outlier_threshold | X < means - sds*outlier_threshold] <- NA

	message("Estimating means")
	means <- rowMeans(X, na.rm=T)

	message("Estimating SDs")
	sds <- rowVars(X, na.rm=T)

	message("Estimating medians")
	medians <- rowMedians(X, na.rm=T)

	message("Counting outliers")
	outliers <- apply(X, 1, function(x) sum(is.na(x)))

	dat <- data.frame(cpg=rownames(X), mean=means, median=medians, sd=sds, outlier=outliers)
	return(dat)
}

message("Generating summary stats of methylation")

meth_summary <- summariseMeth(norm.beta, 5)

save(cohort_summary, file=cohort_descriptives_file)
save(meth_summary, file=methylation_summary_file)

message("You successfully performed all datachecks!")
