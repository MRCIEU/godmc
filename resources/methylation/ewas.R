library(matrixStats)

# Perform EWAS


main <- function()
{

	arguments <- commandArgs(T)

	beta_file <- arguments[1]
	phen_file <- arguments[2]
	out_file <- arguments[3]

	load(beta_file)
	phen <- read.table(phen_file, he=T)
	rownames(phen) <- phen$IID
	phen <- subset(phen, IID %in% colnames(norm.beta), select=-c(IID))

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(phen)]
	phen <- phen[match(colnames(norm.beta), rownames(phen)), ]
	stopifnot(all(rownames(phen) == colnames(norm.beta)))

	norm.beta <- t(norm.beta)

	for(i in 1:ncol(phen))
	{
		phen_name <- colnames(phen)[i]
		message("Performing EWAS for ", phen_name)
		res <- perform_assoc(norm.beta, phen[,i])
		message("Saving results")
		save(res, file=paste0(out_file, ".", phen_name, ".RData"))
	}
}


perform_assoc <- function(X, y)
{
	require(matrixStats)
	if(any(is.na(y)))
	{
		message("Removing missing values")
		index <- is.na(y)
		y <- y[!index]
		X <- X[!index, ]
	}
	message("Calcualting correlations")
	r <- cor(X, y)
	message("Generating betas")
	sdy <- sd(y, na.rm=T)
	sdx <- colSds(X)
	b <- r * sdy / sdx

	message("Generating test statistic")
	# n <- apply(X, 2, function(x) sum(!is.na(x)))
	n <- nrow(X) - colSums(is.na(X))
	fstat <- r^2 / ((1 - r^2) / (n - 2))
	pval <- pf(fstat, 1, n-1, low=F)
	tval <- qt(pval/2, n-1, low=FALSE)

	message("Estimating standard errors")
	se <- abs(b) / abs(tval)
	dat <- data.frame(b=b, se, n=n, pval=pval)
	return(dat)
}

main()