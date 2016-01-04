errorlist <- list()
warninglist <- list()


library(data.table)
suppressMessages(library(matrixStats))

args <- commandArgs(TRUE)
cnv_file <- as.character(args[1])
meth_ids_file <- as.character(args[2])
covariates_file <- as.character(args[3])
cohort_descriptives_file <- as.character(args[4])



message("Checking CNV data: ", cnv_file)
load(cnv_file)

if(! "cnv" %in% ls())
{
	msg <- paste0("please save cnv matrix with the object name 'cnv'")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

# is cnv data a matrix

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

meth_ids <- scan(meth_ids_file, what="character")
covar <- read.table(covariates_file, header=TRUE, stringsAsFactors=FALSE)

commonids_mcc <- Reduce(intersect, list(colnames(cnv), covar$IID, meth_ids))
message(length(commonids_mcc), " samples in common between CNV, covariate and methylation data")

if(length(commonids_mcc) < 50)
{
	msg <- paste0("fewer than 50 subjects with methylation, CNV and covariate data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

if(d2 != length(meth_ids))
{
	msg <- paste0("number of subject with CNV data is not the same as number of subjects with methylation data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}


cohort_summary <- list()
cohort_summary$mcnv_sample_size <- length(commonids_mcc)
cohort_summary$n_cnv_positions <- nrow(cnv)


save(cohort_summary, file=cohort_descriptives_file)

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