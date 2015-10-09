main <- function()
{
	arguments <- commandArgs(T)

	cellcount_file <- arguments[1]
	fam_file <- arguments[2]
	out_file <- arguments[3]

	cellcounts <- read.table(cellcount_file, he=T)
	fam <- read.table(fam_file)[,1:2]

	cellcounts$entropy <- apply(as.matrix(cellcounts[,-1]), 1, function(x)
	{
		h <- x * log2(x)
		h[!is.finite(h)] <- 0
		-sum(h)
	})

	for(i in 2:ncol(cellcounts))
	{
		cellcounts[,i] <- rntransform(cellcounts[,i])
	}
	nom <- names(cellcounts)[-1]
	cellcounts <- merge(cellcounts, fam, by.x="IID", by.y="V2")
	cellcounts <- subset(cellcounts, select=c(V1, V2, nom))

	write.table(cellcounts, file=out_file, row=F, col=F, qu=F)	
}



rntransform <- function(x)
{
	out <- rank(x) - 0.5
	out[is.na(x)] <- NA
	out <- out/(max(out, na.rm = T) + 0.5)
	out <- scale(qnorm(out))
	out
}


main()
