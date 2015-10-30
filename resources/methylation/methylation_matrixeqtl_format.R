arguments <- commandArgs(T)

methylationfile <- arguments[1]
sqfile <- arguments[2]

message("Loading methylation data")
load(paste0(methylationfile, ".RData"))

message("Writing MatrixEQTL format")
write.table(norm.beta, file=paste0(methylationfile, ".txt"), row=T, col=T, qu=F, sep="\t")

message("Writing squared z-scores")
norm.beta <- norm.beta^2
write.table(norm.beta, file=paste0(sqfile, ".txt"), row=T, col=T, qu=F, sep="\t")

message("Successfully converted methylation data!")