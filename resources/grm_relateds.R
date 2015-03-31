#!/usr/bin/R

#' Read binary GRM files into R
#'
#' @param rootname
#' @export
#' @return List of GRM and id data frames
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
 
#' Write readGRM style output back to binary GRM for use with GCTA
#'
#' @param grm Output from \link{readGRM}
#' @param rootname
#' @export
writeGRM <- function(grm, rootname)
{
	bin.file.name <- paste(rootname, ".grm.bin", sep="")
	n.file.name <- paste(rootname, ".grm.N.bin", sep="")
	id.file.name <- paste(rootname, ".grm.id", sep="")
	write.table(grm$id, id.file.name, row=F, col=F, qu=F)
	n <- dim(grm$id)[1]
	bin.file <- file(bin.file.name, "wb")
	writeBin(grm$grm$grm, bin.file, size=4)
	close(bin.file)
	n.file <- file(n.file.name, "wb")
	writeBin(grm$grm$N, n.file, size=4)
	close(n.file)
}
 
 
setUnrelsZero <- function(grm, threshold)
{
	index <- grm$grm$grm < threshold
	grm$grm$grm[index] <- 0
	return(grm)
}
 
 
arguments <- commandArgs(T)
 
infile <- arguments[1]
outfile <- arguments[2]
threshold <- as.numeric(arguments[3])
 
grm <- readGRM(infile)
grm <- setUnrelsZero(grm, threshold)
writeGRM(grm, outfile)