library(ggplot2)
library(data.table)

main <- function()
{
	arguments <- commandArgs(T)

	in_file <- arguments[1]
	pval_column <- as.numeric(arguments[2])
	chr_column <- as.numeric(arguments[3])
	pos_column <- as.numeric(arguments[4])
	header <- as.logical(arguments[5])
	control_chr <- as.numeric(arguments[6])
	control_pos <- as.numeric(arguments[7])
	control_window <- as.numeric(arguments[8])
	control_threshold <- as.numeric(arguments[9])
	out <- arguments[10]

	message("Reading in GWAS results")
	a <- as.data.frame(fread(paste("zcat", in_file), header=header), stringsAsFactors=FALSE)

	if(length(unique(a[,chr_column])) > 30)
	{
		stop("Wrong chromosome column specified")
	}

	if(min(a[,pval_column], na.rm=TRUE) < 0 | max(a[,pval_column], na.rm=TRUE) > 1)
	{
		stop("Wrong column specified for p-values")
	}

	if(any(a[, pos_column] < 0))
	{
		stop("Negative values in position column")
	}

	message("Generating QQ-plot")
	lambda <- qqplot_pval(a[,pval_column], plot=TRUE, filename=paste0(out, "_qqplot.png"))

	message("Generating Manhattan plot")
	manhattan_plot(a[,pval_column], a[,chr_column], a[,pos_column], filename=paste0(out, "_manhattan.png"))

	message("The following plots have been generated, please check!\n",
		paste0(out, "_manhattan.png\n"),
		paste0(out, "_qqplot.png")
	)

	index <- a[,pos_column] > (control_pos - control_window) & a[,pos_column] < (control_pos + control_window)

	a <- a[index, ]
	min_pval <- min(a[,pval_column], na.rm=TRUE)

	message("\n\nExpecting a large meQTL near ", control_chr, ":", control_pos)
	message("Lowest p-value within ", control_window, " base pairs:")
	message(min_pval)
	if(min_pval > control_threshold)
	{
		message("WARNING!")
		message("There doesn't appear to be a QTL for this positive control")
		message("Please upload this section and contact GoDMC analysts before continuing.\n\n")
	}

	message("\n\nlambda value for GWAS: ", lambda$estimate)
	if(lambda$estimate > 1.1)
	{
		message("WARNING!")
		message("The median lambda value is higher than expected here")
		message("This suggests population stratification is affecting the GWAS")
		message("Please upload this section and contact GoDMC analysts before continuing.\n\n")
	}

}

manhattan_plot <- function(p, chr, pos, filename=NULL, width=15, height=7, threshold=-log10(0.05/1000000), maxval=20)
{
	dat <- data.frame(chrom=as.numeric(chr), bp=pos, pval=pmin(maxval, -log10(p)))
	dat <- dat[order(dat$chrom, dat$bp), ]
	dat$col <- dat$chr %% 2 + 1
	dat <- subset(dat, !is.na(pval))

	pl <- ggplot(dat, aes(x=bp, y=pval)) +
	geom_point(aes(colour=factor(col))) +
	facet_grid(. ~ chrom, scale="free_x", space="free_x") +
	theme(legend.position="none") +
	scale_colour_manual(values=c("#404040", "#ca0020")) +
	theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
	ylim(0, max(c(threshold, dat$pval, na.rm=TRUE))) +
	labs(y=expression(-log[10]*p), x="Position") +
	geom_hline(yintercept=threshold)

	if(!is.null(filename))
	{
		ggsave(filename, pl, width=width, height=height)
	} else {
		print(pl)
	}
}

qqplot_pval <- function (data, plot = FALSE, filename=NULL, proportion = 1, method = "regression", filter = TRUE, df = 1, ...)
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
	if(plot)
	{
		if(!is.null(filename))
		{
			png(filename)
		}
		lim <- c(0, max(data, ppoi, na.rm = TRUE))
		oldmargins <- par()$mar
		par(mar = oldmargins + 0.2)
		plot(ppoi, data, xlab = expression("Expected " ~ chi^2), 
			ylab = expression("Observed " ~ chi^2), 
			main=paste0(
				"lambda = ", round(out$estimate, 3)),
			...
		)
		abline(a = 0, b = 1)
		abline(a = 0, b = out$estimate, col = "red")
		par(mar = oldmargins)
		if(!is.null(filename))
		{
			dev.off()
		}
	}
	out
}

main()
