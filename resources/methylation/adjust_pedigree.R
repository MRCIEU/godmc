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
	covs <- read.table(cov_file, he=T,stringsAsFactors=F)
	index <- apply(covs, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
	covs <- covs[!index, ]
	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))
	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]

	g <- grep("_factor",names(covs))
	if(length(g) > 0)
	{
		for (i in 1:length(g))
		{
			covs[,g[i]]<-as.factor(covs[,g[i]])
		}
	}
	grm <- readGRM(grmfile)
	kin <- makeGRMmatrix(grm)
	kin <- kin[rownames(kin) %in% colnames(norm.beta), colnames(kin) %in% colnames(norm.beta)]

	message("Kinship matrix is ", nrow(kin), " by ", nrow(kin))

	index <- match(rownames(kin), colnames(norm.beta))
	norm.beta <- norm.beta[,index]
	covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
	stopifnot(all(rownames(kin) == colnames(norm.beta)))
	stopifnot(all(rownames(covs) == colnames(norm.beta)))

	message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	norm.beta.copy <- norm.beta
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta.copy, na.rm=T)
		means <- rowMeans(norm.beta.copy, na.rm=T)
		norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	}
	outlier_count <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	norm.beta.copy <- is.na(norm.beta.copy)

	message("Calculating eigenvectors")

	relmat <- kin * 2
	tmp <- t(relmat)
	relmat[upper.tri(relmat)] <- tmp[upper.tri(tmp)]
	eig <- eigen(relmat, symmetric=TRUE)
	message(class(eig))
	print(str(eig))

	rm(tmp, relmat)


	message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

	if(is.na(nthreads) | nthreads == 1)
	{
		out <- adjust.relatedness.serial(norm.beta, covs, kin, eig)
	} else {
		message("Running with ", nthreads, " threads")
		out <- adjust.relatedness(norm.beta, covs, kin, eig, nthreads)
	}

	norm.beta <- out$x
	classes <- data.frame(cpg=rownames(norm.beta), cl=out$cl)
	norm.beta[norm.beta.copy] <- NA

	index <- which(is.na(norm.beta), arr.ind = TRUE) 

	if (length(index)>0){
    message("Replace ",length(index)," missing values with rowmeans")
    norm.beta[index] <- rowMeans(norm.beta, na.rm = TRUE)[index[, "row"]] }

	save(norm.beta, file=out_file)
	save(classes, file=paste0(out_file, "_classes"))
	message("Successfully completed adjustments")
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


adjust.relatedness.fast.1 <- function(x, covs, kin, eig, quiet=TRUE)
{
	x[!is.finite(x)] <- mean(x, na.rm=T)
	d <- data.frame(X=rntransform(x), covs)
	rownames(d) <- colnames(kin)
	form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
	d$X <- residuals(lm(form, d))
	p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=quiet))

	iter <- 1
	while(class(p_out) == "try-error" & iter < 5)
	{
		message("trying again...")
		iter <- iter + 1
		p_out <- try(polygenic(X, data=d, kinship.matrix=kin, eigenOfRel=eig, quiet=quiet))
	}
	if(class(p_out) == "try-error")
	{
		message("giving up, just using fixed effects model")
		a <- d$X
	} else {
		a <- as.numeric(rntransform(p_out$grresidualY))
	}
	return(list(x=a, cl=class(p_out)))
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

adjust.relatedness <- function(B, covs, kin, eig, mc.cores=mc.cores)
{

	l1 <- get.index.list(nrow(B), mc.cores)
	l <- lapply(l1, function(ii)
	{
		res <- mclapply(ii, function(i)
		{
			# if( i %% 100 == 0) message("Probe ", i, " of ", nrow(B))
			message("Probe ", i, " (", rownames(B)[i], ") of ", nrow(B))
			out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig)
		}, mc.cores=mc.cores)
		a <- do.call(rbind, lapply(res, function(x) x$x))
		b <- sapply(res, function(x) x$cl)
		return(list(x=a, cl=b))
	})
	x <- do.call(rbind, lapply(l, function(x) x$x))
	cl <- unlist(lapply(l, function(x) x$cl))
	rownames(x) <- rownames(B)
	colnames(x) <- colnames(B)
	return(list(x=x, cl=cl))
}



adjust.relatedness.fast <- function(B, covs, kin, eig, mc.cores=mc.cores)
{

	res <- mclapply.safe(1:nrow(B), function(i)
	{
		# if( i %% 100 == 0) message("Probe ", i, " of ", nrow(B))
		message("Probe ", i, " (", rownames(B)[i], ") of ", nrow(B))
		adjust.relatedness.fast.1(B[i,], covs, kin, eig)
	}, mc.cores=mc.cores)

	res <- do.call(rbind, res)
	rownames(res) <- rownames(B)
	colnames(res) <- colnames(B)
	return(res)
}




adjust.relatedness.serial <- function(B, covs, kin, eig)
{
	cl <- array(0, nrow(B))
	for(i in 1:nrow(B))
	{
		cat(i, "\n")
		out <- adjust.relatedness.fast.1(B[i,], covs, kin, eig)
		B[i, ] <- out$x
		cl[i] <- out$cl
	}
	return(list(x=B, cl=cl))

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

