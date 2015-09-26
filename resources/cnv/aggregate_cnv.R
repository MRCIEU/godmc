arguments <- commandArgs(T)
fileroot <- arguments[1]

files <- list.files(pattern=paste(fileroot, "*", sep=""))

message("Aggregating the following", length(files), "files")
message(paste(files, collapse="\n"))

cnv <- list()
for(i in 1:length(files))
{
	message("Reading", files[i])
	load(files[i])
	cnv <- c(cnv, cnv_partial)
}

save(cnv, file=paste(fileroot, ".RData", sep=""))
