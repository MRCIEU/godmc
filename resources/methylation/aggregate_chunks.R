arguments <- commandArgs(T)

rootname <- arguments[1]
n_chunks <- as.numeric(arguments[2])

message("Checking all files are present")
fs <- paste0(rootname, ".", 1:n_chunks, ".RData")
index <- file.exists(fs)
if(any(!index))
{
	message("The following files don't exist, please generate them before continuing.\n", paste(fs[!index], collapse="\n"))
	q()
}
load(fs[1])

message("Reading in data")
alldat <- matrix(NA, nrow(norm.beta) * n_chunks, ncol(norm.beta))
alldat[1:nrow(norm.beta),] <- norm.beta
colnames(alldat) <- colnames(norm.beta)
cpg <- list()
cpg[[1]] <- rownames(norm.beta)
j <- nrow(norm.beta)
for(i in 2:n_chunks)
{
	load(fs[i])
	message(i, " of ", n_chunks)
	i1 <- j + 1
	i2 <- j + nrow(norm.beta)
	alldat[i1:i2, ] <- norm.beta
	cpg[[i]] <- rownames(norm.beta)
	j <- i2
}

norm.beta <- alldat[1:j,]
rownames(norm.beta) <- unlist(cpg)
save(norm.beta, file=paste0(rootname, ".RData"))
message("Successfully aggregated all chunks")
