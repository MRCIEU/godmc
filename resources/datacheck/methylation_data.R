errorlist <- list()
warninglist <- list()


library(data.table)
suppressMessages(library(matrixStats))

args <- (commandArgs(TRUE));
betas_file <- as.character(args[1]);
fam_file <- as.character(args[2]);
cellcounts_file <- as.character(args[3]);
meth_ids_file <- as.character(args[4]);
cohort_descriptives_file <- as.character(args[5])
methylation_summary_file <- as.character(args[6])
ids <- as.character(args[7]);
ids_plink <- as.character(args[8]);


message("Checking methylation data: ", betas_file)
load(betas_file)


fam <- read.table(fam_file, header=FALSE, stringsAsFactors=FALSE)

if(! "norm.beta" %in% ls())
{
	msg <- paste0("please save methylation matrix with the object name norm.beta")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("norm.beta object found")
# is methylation data a matrix

if(!is.matrix(norm.beta))
{
	msg <- paste0("please transform methylation norm.beta to a matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
# are CpGs in rows?
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

# are individuals unique
if(any(duplicated(colnames(norm.beta))))
{
	msg <- paste0("please remove duplicate samples from methylation data")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("No duplicate IDs")

# check for NAs in beta matrix
if(any(is.na(norm.beta)))
{
	msg <- paste0("please remove NAs from methylation matrix")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}
message("No NAs in data")

# check for negative values in beta matrix
if(any(norm.beta < 0))
{
	msg <- paste0("please remove negative values from methylation matrix. Are these beta values?")
	errorlist <- c(errorlist, msg)
	warning("ERROR: ", msg)
}

# check for values above 1 in beta matrix
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

# extract list of individuals with geno+methylation data
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


# CELLCOUNTS

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


cohort_summary <- list()
cohort_summary$n_CpGs <- nrow(norm.beta)
cohort_summary$methylation_sample_size <- ncol(norm.beta)
cohort_summary$geno_meth_common_ids <- length(overlap)

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
write.table(colnames(norm.beta), file=meth_ids_file, row=F, col=F, qu=F)


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
