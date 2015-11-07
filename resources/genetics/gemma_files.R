
main <- function()
{

	arguments <- commandArgs(T)
	grm_file <- arguments[1]
	cellcounts_file <- arguments[2]
	cellcounts_annotation <- arguments[3]
	smoking_file <- arguments[4]

	message("Converting GRM into GEMMA format")
	grm <- readGRM(grm_file)
	ids <- grm$id
	grm <- makeGRMmatrix(grm)

	message("Writing GRM")
	write.table(round(grm, 5), file=paste0(grm_file, ".gemma"), row=F, col=F, qu=F, sep="\t")

	message("Adjusting cellcounts for smoking and converting into GEMMA format")
	cellcounts <- read.table(cellcounts_file, he=T, stringsAsFactors=FALSE)
	cellcounts <- cellcounts[match(cellcounts$IID, ids$V2), ]

	smok <- read.table(smoking_file, header=TRUE)
	smok <- smok[match(smok$IID, cellcounts$IID), ]
	stopifnot(all(smok$IID == cellcounts$IID))

	entropy <- cellcounts$entropy
	cellcounts <- subset(cellcounts, select=-c(IID, entropy))

	# Adjust for smoking
	for(i in 1:ncol(cellcounts))
	{
		cellcounts[,i] <- residuals(lm(cellcounts[,i] ~ smok$Smoking, na.action=na.exclude))
	}

	# Find variables with large number of missing values
	#index <- apply(cellcounts, 2, function(x)
	#{
#		(sum(is.na(x)) / length(x)) < 0.2
	#})

	#cellcounts <- cellcounts[, index]
	#write.table(data.frame(ids, cellcounts), file=paste0(cellcounts_file, ".plink"), row=F, col=F, qu=F)

	# Replace missing values with mean values
	#for(i in 1:ncol(cellcounts))
	#{
#		cellcounts[is.na(cellcounts[,i]), i] <- mean(cellcounts[,i], na.rm=T)
	#}

	write.table(round(cellcounts, 5), file=paste0(cellcounts_file, ".gemma"), row=F, col=F, qu=F, sep="\t")
	write.table(data.frame(ids, entropy), file=paste0(cellcounts_file, ".entropy.plink"), row=F, col=F, qu=F)
	write.table(names(cellcounts), file=cellcounts_annotation, row=F, col=F, qu=F)

}


makeGRMmatrix <- function(grm)
{
	mat <- diag(nrow(grm$id))
	mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
	mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
	return(mat)
}


readGRM <- function(rootname)
{
	bin.file.name <- paste(rootname, ".grm.bin", sep="")
	n.file.name <- paste(rootname, ".grm.N.bin", sep="")
	id.file.name <- paste(rootname, ".grm.id", sep="")
 
	cat("Reading IDs\n")
	id <- read.table(id.file.name)
	n <- dim(id)[1]
	cat("Reading GRM\n")
	bin.file <- file(bin.file.name, "rb")
	grm <- readBin(bin.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(bin.file)
	cat("Reading N\n")
	n.file <- file(n.file.name, "rb")
	N <- readBin(n.file, n=n*(n+1)/2, what=numeric(0), size=4)
	close(n.file)
 
	cat("Creating data frame\n")
	l <- list()
	for(i in 1:n)
	{
		l[[i]] <- 1:i
	}
	col1 <- rep(1:n, 1:n)
	col2 <- unlist(l)
	grm <- data.frame(id1=col1, id2=col2, N=N, grm=grm)	
 
	ret <- list()
	ret$grm <- grm
	ret$id <- id
	return(ret)
}

main()
