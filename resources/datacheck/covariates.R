errorlist <- list()
warninglist <- list()

library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
covariates_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
meth_ids_file <- as.character(args[3])
age_distribution_plot <- as.character(args[4])
cohort_descriptives_file <- as.character(args[5])


message("Checking covariates file: ", covariates_file)
covar <- read.table(covariates_file,header=T)
cov1 <- dim(covar)[1]
cov2 <- dim(covar)[2]

meth_ids <- scan(meth_ids_file, what="character")
fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

commonids_mgc <- Reduce(intersect, list(meth_ids, covar$IID, fam[,2]))
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
if(length(which(a>0.1*length(meth_ids))))
{
	msg <- paste0("more than 10% of missingness in the following covariates:\n", 
		paste(names(covar)[a>0.1*length(meth_ids)], collapse="\n")
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


covar <- subset(covar, IID %in% commonids_mgc)

cohort_summary <- list()
cohort_summary$mqtl_sample_size <- length(commonids_mgc)
cohort_summary$mqtl_n_males <- sum(covar$Sex_factor == "M",na.rm=T)
cohort_summary$mqtl_n_females <- sum(covar$Sex_factor == "F",na.rm=T)
cohort_summary$mqtl_mean_age <- mean(covar$Age_numeric,na.rm=T)
cohort_summary$mqtl_median_age <- median(covar$Age_numeric,na.rm=T)
cohort_summary$mqtl_sd_age <- sd(covar$Age_numeric,na.rm=T)
cohort_summary$mqtl_max_age <- max(covar$Age_numeric,na.rm=T)
cohort_summary$mqtl_min_age <- min(covar$Age_numeric,na.rm=T)
cohort_summary$covariates <- names(covar)[-1]


save(cohort_summary, file=cohort_descriptives_file)


message("\n\nCompleted checks\n")

message("Summary of data:\n")
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
