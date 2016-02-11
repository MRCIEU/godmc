suppressMessages(library(matrixStats))
suppressMessages(library(meffil))

# Perform EWAS


main <- function()
{

	arguments <- commandArgs(T)

	beta_file <- arguments[1]
	phen_file <- arguments[2]
	covs_file <- arguments[3]
	phen_name <- arguments[4]
	min_age <- as.numeric(arguments[5])
	max_age <- as.numeric(arguments[6])
	pc_file <-arguments[7]
	out_file <- arguments[8]
	qqplot_file <- arguments[9]
    

	message("Loading methylation data")
	phen <- read.table(phen_file, he=T, stringsAsFactors=FALSE)
	if(! phen_name %in% names(phen))
	{
		message(phen_name, " not available. Stopping the analysis.")
		q()
	}
    rownames(phen) <- phen$IID

    phen<-phen[,which(names(phen)%in%c("IID",phen_name))]

	covs <- read.table(covs_file, he=T, stringsAsFactors=FALSE)
    g<-grep("factor",names(covs))
    if(length(g)>1){ 
    for (i in 1:length(g)){
    covs[,g[i]]<-as.factor(covs[,g[i]])
    if(length(levels(covs[,g[i]]))==1)
    covs<-covs[,-g[i]]
    }
    }

    g<-grep("numeric",names(covs))
    if(length(g)>1){ 
    for (i in 1:length(g)){
    covs[,g[i]]<-as.numeric(as.character(covs[,g[i]]))
    }
    }

	rownames(covs) <- covs$IID
	
    pcs<-read.table(paste(pc_file,".txt",sep=""),sep=" ",header=T)
    rownames(pcs) <- pcs$IID
    m<-match(covs$IID,pcs$IID)
    covs<-data.frame(covs,pcs[m,-1])
  
	if(nrow(covs) < 10)
	{
		message("There are fewer than 10 individuals remaining. Stopping the analysis.")
		q()
	}

	load(beta_file)

    keep_ids <- subset(covs, Age_numeric >= min_age & Age_numeric < max_age)$IID

	phen <- subset(phen, IID %in% colnames(norm.beta) & IID %in% keep_ids, select=phen_name)
	if(sum(!is.na(phen[[phen_name]])) < 10)
	{
		message("There are fewer than 10 individuals remaining. Stopping the analysis.")
		q()
	}
    
    phen<-na.omit(phen)

	norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(phen)]
	phen <- phen[match(colnames(norm.beta), rownames(phen)), , drop=FALSE]
	stopifnot(all(rownames(phen) == colnames(norm.beta)))
    
    covs <- covs[match(colnames(norm.beta), rownames(covs)), , drop=FALSE]
	covs<-covs[,which(colnames(covs)%in%c("IID")==F)]
    
    message("\nPerforming EWAS for ", phen_name)
    ewas.ret <- meffil.ewas(norm.beta, variable=phen[,1], covariates=covs,winsorize.pct = NA,most.variable = min(nrow(norm.beta), 20000))

    res<-as.data.frame(ewas.ret$p.value)
  
	message("Generating Q-Q plot")
	qqplot_pval(res$none, file=paste(qqplot_file,"nocovs.png",sep="."))
    qqplot_pval(res$all, file=paste(qqplot_file,"allcovs.png",sep="."))
    qqplot_pval(res$isva0, file=paste(qqplot_file,"isva0.png",sep="."))
    qqplot_pval(res$isva1, file=paste(qqplot_file,"isva1.png",sep="."))	

	message("Saving results")
	save(ewas.ret, file=out_file)

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