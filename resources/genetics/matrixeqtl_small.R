library(dplyr)
library(data.table)
library(MatrixEQTL)


main <- function()
{
	arguments <- commandArgs(T)
	geno_file <- arguments[1]
	phen_file <- arguments[2]
	cov_file <- arguments[3]
	threshold <- as.numeric(arguments[4])
	snplist <- arguments[5]
	out_file <- arguments[6]

	snplist <- fread(paste0("zcat ", snplist), header=TRUE)
	snplist <- subset(snplist, select=c(SNP, A1, A2, MAF))
	names(snplist) <- c("id", "EA", "NEA", "EAF")
	snplist$num <- 1:nrow(snplist)

	slicesize <- 1000

	useModel = modelLINEAR
	errorCovariance = numeric()

	snps = SlicedData$new()
	snps$fileDelimiter = "\t"
	snps$fileOmitCharacters = "NA"
	snps$fileSkipRows = 1
	snps$fileSkipColumns = 1
	snps$fileSliceSize = slicesize
	snps$LoadFile( geno_file )

	message("Number of SNPs: ", nrow(snplist), " ", nrow(snps))

	stopifnot(nrow(snplist) == nrow(snps))
	stopifnot(all(rownames(snps) == snplist$id))

	cvrt <- SlicedData$new()
	if(cov_file!="NULL")
	{
		cvrt$fileDelimiter = "\t"
		cvrt$fileOmitCharacters = "NA"
		cvrt$fileSkipRows = 1
		cvrt$fileSkipColumns = 1
		cvrt$fileSliceSize = slicesize
		cvrt$LoadFile( cov_file )
	}

	gene <- SlicedData$new()
	gene$fileDelimiter = "\t"
	gene$fileOmitCharacters = "NA"
	gene$fileSkipRows = 1
	gene$fileSkipColumns = 1
	gene$fileSliceSize = slicesize
	gene$LoadFile( phen_file )

	cpglist <- data.frame(id=rownames(gene), num=1:nrow(gene), stringsAsFactors=FALSE)

	if(nrow(cvrt) > 0)
	{
		ids <- Reduce(intersect, list(snps$columnNames, gene$columnNames, cvrt$columnNames))
		snps$ColumnSubsample(match(ids, snps$columnNames))
		gene$ColumnSubsample(match(ids, gene$columnNames))
		cvrt$ColumnSubsample(match(ids, cvrt$columnNames))
		stopifnot(all(snps$columnNames==gene$columnNames))
		stopifnot(all(snps$columnNames==cvrt$columnNames))
	} else {
		ids <- Reduce(intersect, list(snps$columnNames, gene$columnNames))
		snps$ColumnSubsample(match(ids, snps$columnNames))
		gene$ColumnSubsample(match(ids, gene$columnNames))
		stopifnot(all(snps$columnNames==gene$columnNames))
	}

	sample_size <- length(intersect(colnames(gene), colnames(snps)))
	message("Performing matrixeqtl analysis on ", sample_size, " samples")

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = NULL,
		pvOutputThreshold = threshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = FALSE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE
	)

	save_results(me$all$eqtls, cpglist, snplist, sample_size, out_file)
}

save_results <- function(dat, cpglist, snplist, sample_size, out_file)
{
	message("Number of associations: ", nrow(dat))
	message("Calculating SE")
	dat$se <- abs(dat$beta / dat$statistic)
	dat <- subset(dat, select=c(snps, gene, beta, se))

	message("Creating expected list")
	d1 <- expand.grid(cpgid=1:nrow(cpglist), snpid=1:nrow(snplist))
	d1$nom <- paste0(d1$cpgid, "_", d1$snpid)

	message("Checking coverage")
	dat$cpgid <- cpglist$num[match(dat$gene, cpglist$id)]
	dat$snpid <- snplist$num[match(dat$snps, snplist$id)]
	dat$nom <- paste0(dat$cpgid, "_", dat$snpid)
	missing <- subset(d1, ! nom %in% dat$nom)

	if(nrow(missing) >= 1)
	{
		message(nrow(missing), " missing associations")
		empty <- data.frame(
			snps=subset(snplist, num %in% missing$snpid)$id,
			gene=subset(cpglist, num %in% missing$cpgid)$id,
			beta=NA,
			se=NA,
			stringsAsFactors=FALSE
		)
		missing <- cbind(empty, missing)
		dat <- bind_rows(dat, missing)
	}

	message("Ordering results")
	dat <- dat[order(dat$snpid, dat$cpgid), ]
	stopifnot(nrow(dat) == nrow(d1))
	stopifnot(all(dat$nom == d1$nom))

	message("Writing")
	unlink(out_file)
	zz <- gzfile(out_file, "wb")

	# write sample size
	writeBin(as.integer(sample_size), con=zz)

	# Write cpg list
	ncpg <- nrow(cpglist)
	writeBin(as.integer(ncpg), con=zz)
	writeBin(as.character(cpglist$id), con=zz)
	
	# write snplist
	nsnp <- nrow(snplist)
	writeBin(as.integer(nsnp), con=zz)
	writeBin(as.character(snplist$id), con=zz)
	writeBin(as.character(snplist$EA), con=zz)
	writeBin(as.character(snplist$NEA), con=zz)
	writeBin(as.numeric(snplist$EAF), con=zz)

	# write beta
	# write se
	writeBin(as.integer(nrow(dat)), con=zz)
	writeBin(as.numeric(dat$beta), con=zz, size=4)
	writeBin(as.numeric(dat$se), con=zz, size=4)

	close(zz)
	message("Done")
}


make_gwama <- function(out_file)
{
	message("Reading from ", out_file)
	zz <- gzfile(out_file, "rb")
	sample_size <- readBin(zz, integer(), 1)
	message("Sample size: ", sample_size)
	ncpg <- readBin(zz, integer(), 1)
	message("Number of CpGs: ", ncpg)
	cpgid <- readBin(zz, character(), ncpg)
	nsnp <- readBin(zz, integer(), 1)
	message("Number of SNPs: ", nsnp)
	snpid <- readBin(zz, character(), nsnp)
	a1 <- readBin(zz, character(), nsnp)
	a2 <- readBin(zz, character(), nsnp)
	maf <- readBin(zz, numeric(), nsnp)
	message("Expected number of associations: ", nsnp * ncpg)
	nres <- readBin(zz, integer(), 1)
	message("Reading number of associations: ", nres)
	beta <- readBin(zz, numeric(), nres, size=4)
	se <- readBin(zz, numeric(), nres, size=4)
	close(zz)

	message("Combining data")
	dat <- expand.grid(cpgid=cpgid, snpid=snpid)
	# dat$MARKERNAME <- paste0(dat$cpgid, "_", dat$snpid)
	dat$EA <- a1[match(dat$snpid, snpid)]
	dat$NEA <- a2[match(dat$snpid, snpid)]
	dat$EAF <- maf[match(dat$snpid, snpid)]
	dat$BETA <- beta
	dat$SE <- se
	dat$N <- sample_size
	# dat <- subset(dat, select=-c(snpid, cpgid))
	return(dat)
}

# gwama <- make_gwama(out_file)
# dat <- me$all$eqtls

# rand <- sample(1:nrow(dat), 1)
# dat$beta[rand]
# subset(gwama, cpgid == dat$gene[rand] & snpid == dat$snps[rand])$BETA

main()

