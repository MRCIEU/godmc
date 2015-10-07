library(parallel)


main <- function()
{
	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	cov_file <- arguments[2]
	out_file <- arguments[3]
	nthreads <- as.numeric(arguments[4])
	chunks <- as.numeric(arguments[5])
	jid <- as.numeric(arguments[6])


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

	message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")

	covs <- read.table(cov_file, he=T)
	index <- apply(covs, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
	covs <- covs[!index, ]
	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]
	covs <- covs[match(colnames(norm.beta), rownames(covs)), ]
	stopifnot(all(rownames(covs) == colnames(norm.beta)))

	if(is.na(nthreads) | nthreads == 1)
	{
		norm.beta <- adjust.covs.serial(norm.beta, covs)
	} else {
		norm.beta <- adjust.covs(norm.beta, covs, nthreads)	
	}
	
	save(norm.beta, file=out_file)
}


adjust.covs.1 <- function(x, covs)
{
	d <- data.frame(X=rntransform(x), covs)
	form <- as.formula(paste0("X ~ ", paste(names(d)[-1], collapse=" + ")))
	rntransform(residuals(lm(form, data=d)))
}

adjust.covs <- function(B, covs, mc.cores=mc.cores)
{
	tmpList = lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(B), by=mc.cores) })

	message("Adjusting data for covariates, may take some time...")	

	tmpAdj <- mclapply(tmpList, function(ix)
	{ 
		apply(B[ix,], 1, function(x) adjust.covs.1(x, covs))
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

adjust.covs.serial <- function(B, covs)
{
	for(i in 1:nrow(B))
	{
		cat(i, "\n")
		B[i, ] <- adjust.covs.1(B[i,], covs)
	}
	return(B)

	# apply(B, 1, function(x) adjust.relatedness.1(x, kin))
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

main()