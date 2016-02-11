library(MatrixEQTL)

main <- function()
{
	arguments <- commandArgs(T)
	geno_file <- arguments[1]
	phen_file <- arguments[2]
	cov_file <- arguments[3]
	threshold <- as.numeric(arguments[4])
	out_file <- arguments[5]

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
	save(me, file=out_file)
}


main()

