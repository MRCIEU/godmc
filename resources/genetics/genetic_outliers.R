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


removeOutliersFromData <- function(tbetas, thresh, iterations)
{
	n <- nrow(tbetas)
	res <- array(0, c(n, iterations+1))
	for(i in 1:n)
	{
		cat(i, "\n")
		res[i,1] <- sum(!is.na(tbetas[i, ]))
		a <- remRec(tbetas[i, ], thresh, iterations)
		res[i, ] <- a$its
		tbetas[i, ] <- a$x
	}
	return(list(tbetas=tbetas, res=res))
}



arguments <- commandArgs(T)
pcafile <- arguments[1]
pcasd <- as.numeric(arguments[2])
npc <- as.numeric(arguments[3])
out_file <- arguments[4]

pca <- read.table(pcafile)

pca2 <- removeOutliersFromData(pca[,-c(1:2)], pcasd, iterations=3)
index <- apply(pca2, 1, function(x) any(is.na(x)))

genetic_outliers <- pca[index,1:2]

write.table(genetic_outliers, file=out_file, row=F, col=F, qu=F)
if(length(genetic_outliers) > 0)
{
	write.table(subset(pca, !V2 %in% genetic_outliers), file=pcafile, row=F, col=F, qu=F)	
}
