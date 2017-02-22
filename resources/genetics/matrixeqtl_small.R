library(data.table)
library(MatrixEQTL)

main <- function()
{
	arguments <- commandArgs(T)
	print(arguments)
	geno_file <- arguments[1]
	phen_file <- arguments[2]
	cov_file <- arguments[3]
	threshold <- as.numeric(arguments[4])
	cpglist <- arguments[5]
	snplist <- arguments[6]
	freq <- arguments[7]
	out_file <- arguments[8]

	freq <- fread(paste0("zcat ", freq))
	cpglist <- scan(cpglist, what="character")
	snplist <- scan(snplist, what="character")
	freq <- subset(freq, V2 %in% snplist, select=c(V2, V3, V4, V5))
	names(freq) <- c("snps", "EA", "NEA", "EAF")
	cpglist <- data.frame(CPG=cpglist, cpgid=1:length(cpglist), stringsAsFactors=FALSE)
	snplist <- data.frame(SNP=snplist, snpid=1:length(snplist), stringsAsFactors=FALSE)


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

	message("Performing matrixeqtl analysis on ", length(snps$columnNames), " samples")

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = NULL,
		pvOutputThreshold = threshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = TRUE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE
	)

	save_results(me, cpglist, snplist, freq, out_file)
}


save_results <- function(dat, cpglist, snplist, freq, out_file)
{
	dat$cpgid <- cpglist$cpgid[match(a$PHENOTYPE, cpglist$CPG)]
	dat$snpid <- snplist$snpid[match(a$PHENOTYPE, snplist$SNP)]
	dat$MARKERNAME <- paste0(dat$cpgid, "_", dat$snpid)
	dat$SE <- abs(dat$beta / dat$statistic)
	dat$BETA <- dat$beta
	dat <- merge(dat, freq, by="snps")
	dat <- subset(dat, select=c(MARKERNAME, EA, NEA, EAF, BETA, SE))
	con <- gzfile(out_file, "w")
	write.table(dat, con, row=FALSE, col=TRUE, qu=FALSE)
}

main()

