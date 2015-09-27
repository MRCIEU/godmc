library(parallel)
library(GenABEL)

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

makeGRMmatrix <- function(grm, ids)
{
	mat <- diag(nrow(grm$id))
	mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
	mat <- t(mat)
	nsnpvec <- subset(grm$grm, id1 != id2)$N
	mat[upper.tri(mat, diag=FALSE)] <- nsnpvec
	rownames(mat) <- grm$id$V2
	colnames(mat) <- grm$id$V2
	mat <- mat[match(ids, rownames(mat)), match(ids, colnames(mat))]
	return(mat)
}

adjust.relatedness.1 <- function(x, kin, quiet=TRUE)
{
	d <- data.frame(x)
	rntransform(polygenic(x, d, kin, quiet=quiet)$pgresidualY)
}

rntransform <- function(x)
{
	out <- rank(x) - 0.5
	out[is.na(x)] <- NA
	mP <- 0.5/max(out, na.rm = T)
	out <- out/(max(out, na.rm = T) + 0.5)
	out <- scale(qnorm(out))
	out
}

adjust.relatedness <- function(B, kin, mc.cores=mc.cores)
{
	tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })

	message("Adjusting data for polygenic effects, may take some time minutes...")
	tmpAdj = mclapply(tmpList, function(ix){ apply(B[ix,], 1, function(x) adjust.relatedness.1(x, kin)) }, mc.cores=mc.cores)

	message("Reducing results...")
	adjBeta = matrix(NA, nrow(B), ncol(B))
	for (i in 1:length(tmpList)){
			adjBeta[tmpList[[i]],] = t(tmpAdj[[i]])
	}
	return(adjBeta)
}





arguments <- commandArgs(T)

ccrnmethdatafile <- arguments[1]
grmfile <- arguments[2]
ccrnfammethdatafile <- arguments[3]
nthreads <- as.numeric(arguments[4])

mbeta <- read.table(ccrnfammethdatafile, he=T)

ids <- colnames(mbeta)
grm <- readGRM(grmfile)
kin <- makeGRMmatrix(grm, ids)

mbeta <- adjust.relatedness(mbeta, kin, nthreads)

save(mbeta, file=ccrnfammethdatafile)
