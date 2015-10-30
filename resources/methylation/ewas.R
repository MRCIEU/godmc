suppressMessages(library(matrixStats))

# Perform EWAS


main <- function()
{

	arguments <- commandArgs(T)

	beta_file <- arguments[1]
	phen_file <- arguments[2]
	out_file <- arguments[3]
	qqplot_file <- arguments[4]

	message("Loading methylation data")
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
		message("\nPerforming EWAS for ", phen_name)
		res <- perform_assoc(norm.beta, phen[,i])
		message("Generating Q-Q plot")
		qqplot_pval(res$pval, file=paste0(qqplot_file, ".", phen_name, ".png"))
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


est_lambda <- function (data, plot = FALSE, proportion = 1, method = "regression", filter = TRUE, df = 1, ...)
{
	data <- data[which(!is.na(data))]
	if (proportion > 1 || proportion <= 0) 
		stop("proportion argument should be greater then zero and less than or equal to one")
	ntp <- round(proportion * length(data))
	if (ntp < 1) 
		stop("no valid measurements")
	if (ntp == 1) {
		warning(paste("One measurement, lambda = 1 returned"))
		return(list(estimate = 1, se = 999.99))
	}
	if (ntp < 10) 
		warning(paste("number of points is too small:", ntp))
	if (min(data) < 0) 
		stop("data argument has values <0")
	if (max(data) <= 1) {
		data <- qchisq(data, 1, lower.tail = FALSE)
	}
	if (filter) {
		data[which(abs(data) < 1e-08)] <- NA
	}
	data <- sort(data)
	ppoi <- ppoints(data)
	ppoi <- sort(qchisq(ppoi, df = df, lower.tail = FALSE))
	data <- data[1:ntp]
	ppoi <- ppoi[1:ntp]
	out <- list()
	if (method == "regression") {
		s <- summary(lm(data ~ 0 + ppoi))$coeff
		out$estimate <- s[1, 1]
		out$se <- s[1, 2]
	}
	else if (method == "median") {
		out$estimate <- median(data, na.rm = TRUE)/qchisq(0.5, 
			df)
		out$se <- NA
	}
	else {
		stop("'method' should be either 'regression' or 'median'!")
	}
	if (plot) {
		lim <- c(0, max(data, ppoi, na.rm = TRUE))
		oldmargins <- par()$mar
		par(mar = oldmargins + 0.2)
		plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
			ylab = expression("Observed " ~ chi^2), ...)
		abline(a = 0, b = 1)
		abline(a = 0, b = out$estimate, col = "red")
		par(mar = oldmargins)
	}
	out
}

qqplot_pval <- function(P, filename=NULL)
{
	l <- est_lambda(P, method="median")
	nom <- paste("lambda = ", round(l$estimate, 3), sep="")
	if(!is.null(filename))
	{
		png(filename)
	}
	est_lambda(P, method="median", plot=TRUE, main=nom)
	if(!is.null(filename))
	{
		dev.off()
	}
}


main()