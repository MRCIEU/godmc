library(parallel)
suppressMessages(library(matrixStats))
suppressMessages(library(GenABEL))

main <- function()
{
	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	grmfile <- arguments[2]
	cov_file <- arguments[3]
	out_file <- arguments[4]
	nthreads <- as.numeric(arguments[5])
	chunks <- as.numeric(arguments[6])
	jid <- as.numeric(arguments[7])


	message("Reading methylation data...")
	load(methylationfile)

	if(!is.na(jid))
	{
		chunksize <- ceiling(nrow(norm.beta) / chunks)
		i1 <- chunksize * (jid-1) + 1
		i2 <- min(nrow(norm.beta), chunksize * jid)
		norm.beta <- norm.beta[i1:i2,]
		out_file <- paste0(out_file, ".", jid, ".RData")
	} else {
		out_file <- paste0(out_file, ".RData")
	}

	# Remove all IDs that have any NAs in the covariate file
	covs <- read.table(cov_file, he=T)
	index <- apply(covs, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
	covs <- covs[!index, ]
	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))
	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]

	grm <- readGRM(grmfile)
	kin <- makeGRMmatrix(grm)
	kin <- kin[rownames(kin) %in% colnames(norm.beta), colnames(kin) %in% colnames(norm.beta)]
	index <- match(rownames(kin), colnames(norm.beta))
	norm.beta <- norm.beta[,index]
	covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
	stopifnot(all(rownames(kin) == colnames(norm.beta)))
	stopifnot(all(rownames(covs) == colnames(norm.beta)))

	message("Setting methylation outliers to missing")

	niter <- 3
	outlier_threshold <- 10
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta, na.rm=T)
		means <- rowMeans(norm.beta, na.rm=T)
		norm.beta[norm.beta > means + sds*outlier_threshold | norm.beta < means - sds*outlier_threshold] <- NA
	}

	outlier_count <- apply(norm.beta, 1, function(x) sum(is.na(x)))
	norm.beta.orig <- norm.beta
	print(table(outlier_count))

	message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

	if(is.na(nthreads) | nthreads == 1)
	{
		norm.beta <- adjust.relatedness.serial(norm.beta, covs, kin)
	} else {
		message("Running with ", nthreads, " threads")
		norm.beta <- adjust.relatedness(norm.beta, covs, kin, nthreads)
	}

	message("Checking for any issues with parallelisation")
	outlier_count2 <- apply(norm.beta, 1, function(x) sum(is.na(x)))
	problems <- which(outlier_count != outlier_count2)
	message("Number of parallelisation issues: ", length(problems))

	if(length(problems) > 0)
	{
		message("Recalculating problem probes serially")
		norm.beta.problems <- adjust.relatedness.serial(norm.beta[problems, ], covs, kin)
		norm.beta[problems, ] <- norm.beta.problems
	}

	save(norm.beta, file=out_file)
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

adjust.relatedness.1 <- function(x, covs, kin, quiet=TRUE)
{
	# x <- remRec(x, 10, 3)$x
	d <- data.frame(X=rntransform(x), covs)
	rownames(d) <- colnames(kin)
	form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
	as.numeric(rntransform(polygenic(form, data=d, kinship.matrix=kin, quiet=quiet)$grresidualY))
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

adjust.relatedness <- function(B, covs, kin, mc.cores=mc.cores)
{

	l1 <- get.index.list(nrow(B), mc.cores)
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(i)
		{
			if( i %% 100 == 0) message("Probe ", i, " of ", nrow(B))
			adjust.relatedness.1(B[i,], covs, kin)
		})
		return(do.call(rbind, res))
	})
	l <- do.call(rbind, l)
	rownames(l) <- rownames(B)
	colnames(l) <- colnames(B)
	return(l)
}



adjust.relatedness.serial <- function(B, covs, kin)
{
	for(i in 1:nrow(B))
	{
		cat(i, "\n")
		B[i, ] <- adjust.relatedness.1(B[i,], covs, kin)
	}
	return(B)

	# apply(B, 1, function(x) adjust.relatedness.1(x, kin))
}

get.index.list <- function(n, mc.cores)
{
	mc.cores <- ifelse(mc.cores < 1, 1, mc.cores)
	div <- floor(n / mc.cores)
	rem <- n %% mc.cores
	l1 <- lapply(1:div, function(x) (x-1) * mc.cores + 1:mc.cores)
	if(rem != 0) l1[[div+1]] <- l1[[div]][mc.cores] + 1:rem
	return(l1)
}

removeOutliers <- function(x, thresh, remove=FALSE)
{
	m <- mean(x, na.rm=T)
	s <- sd(x, na.rm=T)
	index <- x > m + thresh*s | x < m - thresh*s

	if(remove)
	{
		x <- x[!index]
	} else {
		x[index] <- NA
	}
	return(x)
}


remRec <- function(x, thresh, iterations)
{
	d <- array(0, iterations+1)
	d[1] <- sum(!is.na(x))
	for(i in 1:iterations)
	{
		x <- removeOutliers(x, thresh)
		d[i+1] <- sum(!is.na(x))
	}
	return(list(x=x, its=d))
}



main()

