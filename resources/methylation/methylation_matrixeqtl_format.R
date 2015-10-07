arguments <- commandArgs(T)

methylationfile <- arguments[1]
sqfile <- arguments[2]

load(paste0(methylationfile, ".RData"))
write.table(norm.beta, file=paste0(methylationfile, ".txt"), row=T, col=T, qu=F, sep="\t")

norm.beta <- norm.beta^2
write.table(norm.beta, file=paste0(sqfile, ".txt"), row=T, col=T, qu=F, sep="\t")
