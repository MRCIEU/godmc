library(CopyNumber450k)
library(CopyNumber450kData)
library(plyr)
library(meffil)
library(splitstackshape)


get_sex <- function(rgset)
{
	gmset <- mapToGenome(preprocessRaw(rgset))
	control_sex <- getSex(gmset)
	pd <- pData(rgset)
	cs <- control_sex$predictedSex[match(rownames(control_sex), rownames(pd))]
	pd$Sample_Sex <- cs
	pData(rgset) <- pd
	return(rgset)
}


map_cnvs <- function(segments)
{
	probenames <- meffil.get.sites()
	probeinfo <- subset(meffil.probe.info(), name %in% probenames & target == "M", select=c(name, chr, pos))
	segments$loc.start <- as.numeric(segments$loc.start)
	segments$loc.end <- as.numeric(segments$loc.end)
	segments$seg.mean <- as.numeric(segments$seg.mean)
	segments$num.mark <- as.numeric(segments$num.mark)
	segments$adjusted.pvalue <- as.numeric(segments$adjusted.pvalue)
	gene <- cSplit(indt=segments, splitCols="genes", sep=";", direction="long")
	gene <- gene[order(gene$genes, gene$adjusted.pvalue), ]
	gene <- subset(gene[order(gene$genes, gene$adjusted.pvalue), ], !duplicated(genes), select=c(genes, seg.mean, num.mark, adjusted.pvalue))
	res <- ddply(segments, .(chrom, loc.start), function(x)
	{
		p <- subset(probeinfo, chr == x$chrom[1] & pos >= x$loc.start[1] & pos <= x$loc.end[1])
		p$seg.mean <- x$seg.mean[1]
		p$num.mark <- x$num.mark[1]
		p$adjusted.pvalue <- x$adjusted.pvalue[1]
		return(p)
	})
	p <- subset(probeinfo, select=c(name))
	p <- merge(p, subset(res, select=c(name, seg.mean, num.mark, adjusted.pvalue)), by="name", all.x=TRUE)
	return(list(cpg=p, gene=gene))
}


calc_cnvs <- function(samplesheet, row, RGcontrolSetEx)
{
	RGset <- read.450k.exp(target=samplesheet[row,])
	RGset <- combine(RGcontrolSetEx, RGset)
	RGset <- get_sex(RGset)
	mcds <- CNV450kSet(RGset)
	mcds <- dropSNPprobes(mcds, maf_threshold=0.01)
	mcds.q <- normalize(mcds, "quantile")
	mcds.q.seg <- segmentize(mcds.q)
	segments <- getSegments(mcds.q.seg)[[1]]
	return(map_cnvs(segments))
}


get_cnvs <- function(samplesheet, RGcontrolSetEx, mc.cores)
{
	tmpList <- lapply(1:mc.cores, function(i){ seq(from=i, to=nrow(samplesheet), by=mc.cores) })

	message("Estimating CNVs from IDAT intensities...")
	tmpAdj <- mclapply(tmpList, function(ix)
	{
		l <- list()
		for(i in ix)
		{
			l[[samplesheet$Sample_Name[i]]] <- calc_cnvs(samplesheet, i, RGcontrolSetEx)
		}
		return(l)
	}, mc.cores=mc.cores)
	tmpAdj <- unlist(tmpAdj, recursive=FALSE)
	return(tmpAdj)
}




arguments <- commandArgs(T)
path <- arguments[1]
outfile <- arguments[2]
ncores <- as.numeric(arguments[3])
jid <- as.numeric(arguments[4])
splitsize <- as.numeric(arguments[5])


# Get idat file info
# path <- "~/data/test_meffil"
samplesheet <- meffil.create.samplesheet(path)
samplesheet$Sample_Group <- "grp1"
samplesheet$Slide <- as.numeric(samplesheet$Slide)

# If this is a batch job then jid will be defined
if(!is.na(jid))
{
	i1 <- (jid - 1) * splitsize + 1
	i2 <- min(jid * splitsize, nrow(samplesheet)
	samplesheet <- samplesheet[i1:i2, ]
}

# Get control data
data(RGcontrolSetEx)

# Get CNVs
cnv <- get_cnvs(samplesheet[1:3,], RGcontrolSetEx, ncores)

# Save
if(!is.na(jid))
{
	cnv_partial <- cnv
	save(cnv_partial, file=outfile)
} else {
	save(cnv, file=outfile)
}

