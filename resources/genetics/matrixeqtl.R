library(MatrixEQTL)

main <- function()
{
	arguments <- commandArgs(T)
	geno_file <- arguments[1]
	phen_file <- arguments[2]
	threshold <- as.numeric(arguments[3])
	out_file <- arguments[4]

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

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = character(),
		output_file_name = out_file,
		pvOutputThreshold = threshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = TRUE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE
	)
	save(me, file=paste(out_file, ".RData", sep=""))
}


main()

