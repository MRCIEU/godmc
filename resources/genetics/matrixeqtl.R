library(MatrixEQTL)

main <- function()
{
	arguments <- commandArgs(T)
	geno_file <- arguments[1]
	phen_file <- arguments[2]
	cov_file <- arguments[3]
	threshold <- as.numeric(arguments[4])
	out_file <- arguments[5]

	useModel = modelLINEAR
	errorCovariance = numeric()

	snps = SlicedData$new()
	snps$fileDelimiter = "\t"
	snps$fileOmitCharacters = "NA"
	snps$fileSkipRows = 1
	snps$fileSkipColumns = 1
	snps$fileSliceSize = 2000
	snps$LoadFile( geno_file )

	gene <- SlicedData$new()
	gene$fileDelimiter = "\t"
	gene$fileOmitCharacters = "NA"
	gene$fileSkipRows = 1
	gene$fileSkipColumns = 1
	gene$fileSliceSize = 2000
	gene$LoadFile( phen_file )

	cvrt <- SlicedData$new()
	cvrt$fileDelimiter = "\t"
	cvrt$fileOmitCharacters = "NA"
	cvrt$fileSkipRows = 1
	cvrt$fileSkipColumns = 1
	cvrt$fileSliceSize = 2000
	cvrt$LoadFile( cov_file )

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = cvrt,
		output_file_name = out_file,
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

