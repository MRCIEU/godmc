suppressPackageStartupMessages(library(meffil))

arguments <- commandArgs(T)

methylationfile <- arguments[1]
methylationfile_pcadj <- arguments[2]
sqfile <- arguments[3]
cov_file <- arguments[4]

message("Loading methylation data")
load(paste0(methylationfile, ".RData"))

message("Writing MatrixEQTL format")
write.table(norm.beta, file=paste0(methylationfile, ".txt"), row=T, col=T, qu=F, sep="\t")

message("Loading pc adjusted methylation data")
load(paste0(methylationfile_pcadj, ".RData"))

message("Writing MatrixEQTL format")
write.table(norm.beta, file=paste0(methylationfile_pcadj, ".txt"), row=T, col=T, qu=F, sep="\t")

message("Writing squared z-scores")
norm.beta <- norm.beta^2
write.table(norm.beta, file=paste0(sqfile, ".txt"), row=T, col=T, qu=F, sep="\t")


# Get males and females separately

covdat <- read.table(cov_file, header=TRUE, stringsAsFactors=FALSE)



male_ids <- subset(covdat, Sex_factor == "M")$IID
xsites <- meffil.get.x.sites()

if(sum(colnames(norm.beta) %in% male_ids > 0))
{
	message("Writing X chromosome methylation data for males")
	write.table(norm.beta[rownames(norm.beta) %in% xsites, colnames(norm.beta) %in% male_ids], file=paste0(methylationfile_pcadj, "_males.txt"), row=T, col=T, qu=F, sep="\t")
}


female_ids <- subset(covdat, Sex_factor == "F")$IID

if(sum(colnames(norm.beta) %in% female_ids > 0))
{
	message("Writing X chromosome methylation data for females")
	write.table(norm.beta[rownames(norm.beta) %in% xsites, colnames(norm.beta) %in% female_ids], file=paste0(methylationfile_pcadj, "_females.txt"), row=T, col=T, qu=F, sep="\t")
}
