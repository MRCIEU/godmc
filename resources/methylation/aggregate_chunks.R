arguments <- commandArgs(T)

rootname <- arguments[1]
n_chunks <- as.numeric(arguments[2])

fs <- paste0(rootname, ".", 1:n_chunks, ".RData")
index <- file.exists(fs)
if(any(!index))
{
	message("The following files don't exist, please generate them before continuing.\n", paste(fs[!index], collapse="\n"))
	q()
}
load(fs[1])

alldat <- matrix(NA, nrow(norm.beta) * n_chunks, ncol(norm.beta))
alldat[1:nrow(norm.beta),] <- norm.beta

j <- nrow(norm.beta)
for(i in 2:n_chunks)
{
	load(fs[i])
	message(i, " of ", n_chunks)
	i1 <- j + 1
	i2 <- j + nrow(norm.beta)
	alldat[i1:i2, ] <- norm.beta
	j <- i2
}

norm.beta <- alldat[1:j,]
save(norm.beta, file=paste0(rootname, ".RData"))