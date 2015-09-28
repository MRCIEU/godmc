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

partition.integer.subsequence <- function(start, end, n) {
    stopifnot(start <= end)
    stopifnot(n <= end-start+1)
    partitions <- floor(seq(start,end+1,length.out=n+1))
    cbind(start=head(partitions, n=-1),
          end=tail(partitions, n=-1) - 1)
}

mclapply.safe <- function (X, FUN, ..., max.bytes=2^30-1) {
    stopifnot(length(X) > 0)

    cores <- options()$mc.cores
    if (is.null(cores) || cores == 1) return(mclapply(X, FUN, ...))
    
    first <- mclapply(X[1], FUN, ...)
    ret.bytes <- as.numeric(object.size(first))
    
    if (length(X) == 1) return(first)

    X <- X[-1]
    
    max.bytes <- min(max.bytes, ret.bytes * length(X))
    n.fun <- floor(max.bytes/ret.bytes)
    if (n.fun < 1)
        stop(paste("The max.bytes parameter is too small.  Try setting it > ", ret.bytes, ".", sep=""))
    
    n.mclapply <- ceiling(length(X)/n.fun)

    partitions <- partition.integer.subsequence(1,length(X),n.mclapply)
    c(first, do.call(c, lapply(1:nrow(partitions), function(i) {
        idx <- partitions[i,"start"]:partitions[i,"end"]
        ret <- mclapply(X[idx], FUN, ...)
        if (length(idx) != length(ret) || any(sapply(ret, is.null)))
            stop(paste("The operating system has decided that some forks of mclapply are using too much memory.\n",
                       "Try reducing the max.bytes parameter or the R option 'mc.cores'."))
        ret
    })))
}

makeGRMmatrix <- function(grm)
{
	mat <- diag(nrow(grm$id))
	mat[upper.tri(mat, diag=TRUE)] <- grm$grm$grm
	mat <- t(mat)
	nsnpvec <- subset(grm$grm, id1 != id2)$N
	mat[upper.tri(mat, diag=FALSE)] <- nsnpvec
	rownames(mat) <- grm$id$V2
	colnames(mat) <- grm$id$V2
	return(mat)
}

adjust.relatedness.1 <- function(x, kin, quiet=TRUE)
{
	d <- data.frame(X=x)
	rownames(d) <- colnames(kin)
	rntransform(polygenic(X, data=d, kinship.matrix=kin, quiet=quiet)$pgresidualY)
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

	message("Adjusting data for polygenic effects, may take some time...")	

	tmpAdj <- mclapply.safe(tmpList, function(ix)
	{ 
		apply(B[ix,], 1, function(x) adjust.relatedness.1(x, kin))
	}, mc.cores=mc.cores)

	message("Reducing results...")
	adjBeta = matrix(NA, nrow(B), ncol(B))
	colnames(adjBeta) <- colnames(B)
	rownames(adjBeta) <- rownames(B)
	for (i in 1:length(tmpList)){
		adjBeta[tmpList[[i]],] = t(tmpAdj[[i]])
	}
	return(adjBeta)
}



adjust.relatedness.serial <- function(B, kin)
{
	# for(i in 1:nrow(B))
	# {
	# 	cat(i, "\n")
	# 	B[i, ] <- adjust.relatedness.1(B[i,], kin)
	# }
	# return(B)

	apply(B, 1, function(x) adjust.relatedness.1(x, kin))
}



arguments <- commandArgs(T)

beta_file <- arguments[1]
grmfile <- arguments[2]
ccrnfammethdatafile <- arguments[3]
nthreads <- as.numeric(arguments[4])
options(mc.cores=nthreads)

load(beta_file)
grm <- readGRM(grmfile)
kin <- makeGRMmatrix(grm)
kin <- kin[rownames(kin) %in% colnames(norm.beta), colnames(kin) %in% colnames(norm.beta)]
index <- match(rownames(kin), colnames(norm.beta))
norm.beta <- norm.beta[,index]
stopifnot(all(rownames(kin) == colnames(norm.beta)))

norm.beta <- adjust.relatedness(norm.beta, kin, nthreads)

write.table(round(norm.beta, 3), file=ccrnfammethdatafile, row=TRUE, col=TRUE, qu=FALSE, sep="\t")
