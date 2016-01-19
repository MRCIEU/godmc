library(plyr)

arguments <- commandArgs(T)

rootname <- arguments[1]
n_chunks <- as.numeric(arguments[2])
out <- arguments[3]

message("Checking all files are present")
fs <- paste0(rootname, ".", 1:n_chunks, ".RData_classes")
index <- file.exists(fs)
if(any(!index))
{
	message("The following files don't exist, please generate them before continuing.\n", paste(fs[!index], collapse="\n"))
	q()
}

message("Reading in data")
l <- list()
for(i in 1:n_chunks)
{
	load(fs[i])
	message(i, " of ", n_chunks)
	l[[i]] <- classes
}

classes <- rbind.fill(l)
save(classes, file=out)
message(sum(classes$cl == "try-error"), " out of ", nrow(classes), " probes couldn't be adjusted for polygenic effects")
message("Successfully aggregated all chunks")
