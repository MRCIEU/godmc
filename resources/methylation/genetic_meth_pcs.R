library(MatrixEQTL)


main <- function()
{
	arguments <- commandArgs(T)
	geno_root <- arguments[1]
	phen_file <- arguments[2]
	out_file <- arguments[3]
	threshold <- as.numeric(arguments[4])

	slicesize <- 1000

	gene <- SlicedData$new()
	gene$fileDelimiter = "\t"
	gene$fileOmitCharacters = "NA"
	gene$fileSkipRows = 1
	gene$fileSkipColumns = 1
	gene$fileSliceSize = slicesize
	gene$LoadFile( paste0(phen_file, ".txt") )

	bn <- basename(geno_root)
	dn <- dirname(geno_root)
	geno_file <- list.files(dn, pattern=paste0("^", bn))

	l <- list()

	for(i in 1:length(geno_file))
	{
		message(i, " of ", length(geno_file))
		snps = SlicedData$new()
		snps$fileDelimiter = "\t"
		snps$fileOmitCharacters = "NA"
		snps$fileSkipRows = 1
		snps$fileSkipColumns = 1
		snps$fileSliceSize = slicesize
		snps$LoadFile( file.path(dn, geno_file[i]) )
		l[[i]] <- run_chunk(snps, gene, threshold, slicesize)
	}
	l1 <- unique(unlist(l))
	load(paste0(phen_file, ".RData"))
	pc <- pc[,!colnames(pc) %in% l1]
	save(pc, file=paste0(out_file, ".RData"))
}




run_chunk <- function(snps, gene, threshold, slicesize)
{

	useModel = modelLINEAR
	errorCovariance = numeric()

	ids <- Reduce(intersect, list(snps$columnNames, gene$columnNames))

	snps$ColumnSubsample(match(ids, snps$columnNames))
	gene$ColumnSubsample(match(ids, gene$columnNames))

	stopifnot(all(snps$columnNames==gene$columnNames))

	me <- Matrix_eQTL_engine(
		snps = snps,
		gene = gene,
		cvrt = SlicedData$new(),
		output_file_name = NULL,
		pvOutputThreshold = threshold,
		useModel = useModel, 
		errorCovariance = errorCovariance, 
		verbose = TRUE,
		pvalue.hist = FALSE,
		min.pv.by.genesnp = FALSE,
		noFDRsaveMemory = FALSE
	)

	return(as.character(unique(me$all$eqtls$gene)))
}


main()

