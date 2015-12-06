library(plyr)
library(meffil)
library(splitstackshape)

arguments <- commandArgs(T)
path <- arguments[1]
outfile <- arguments[2]
ncores <- as.numeric(arguments[3])
jid <- as.numeric(arguments[4])
splitsize <- as.numeric(arguments[5])

options(mc.cores=ncores)

# Get idat file info
# path <- "~/data/test_meffil"
samplesheet <- meffil.create.samplesheet(path)

# If this is a batch job then jid will be defined
if(!is.na(jid))
{
	i1 <- (jid - 1) * splitsize + 1
	i2 <- min(jid * splitsize, nrow(samplesheet))
	samplesheet <- samplesheet[i1:i2, ]
}

library(minfi)
library(CopyNumber450kData)
meffil.add.copynumber450k.references()

segment.lists <- meffil.calculate.cnv(samplesheet, cnv.reference="copynumber450k", chip="450k", verbose=TRUE)
cnv <- meffil.cnv.matrix(segment.lists, "450k")

# Save
if(!is.na(jid))
{
	cnv_partial <- cnv
	save(cnv_partial, file=outfile)
} else {
	save(cnv, file=outfile)
}

