suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))

linear_regression <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)
	p <- pf(fval, 1, n-2, lowe=F)
	return(c(
		bhat=bhat, se=se, pval=p, n=n
	))
}

arguments <- commandArgs(T)

snp_file <- arguments[1]
cpg_file <- arguments[2]
out_file <- arguments[3]
threads <- arguments[4]

message("Reading Inversion data")
snp <- read.table(snp_file)

message("Reading CPG data")
cpg <- fread(cpg_file)


# CPG
message("Organising CpG data")

cpgids <- cpg$V1
cpg <- cpg[, V1:=NULL]
cpg <- t(as.matrix(cpg))
colnames(cpg) <- cpgids

# SNP

message("Organising SNP data")

snp <- t(snp)
colnames(snp) <- snp[1, ]
snp <- snp[-1, ]
rownames(snp) <- snp[, 1]
snp <- snp[, -1]
class(snp) <- "numeric"


snp <- snp[rownames(snp) %in% rownames(cpg), ]
cpg <- cpg[match(rownames(snp), rownames(cpg)), ]
stopifnot(all(rownames(cpg) %in% rownames(snp)))
stopifnot(all(rownames(snp) %in% rownames(cpg)))
stopifnot(identical(rownames(snp), rownames(cpg)))

snpinfo <- data.frame(MARKERNAME=colnames(snp), SNP=colnames(snp), EA = "I", NEA = "N",
                      EAF = colMeans(snp, na.rm = TRUE)/2,
                      stringsAsFactors=FALSE)



message("\n\n")
message("Number of individuals: ", nrow(cpg))
message("Number of CpGs: ", ncol(cpg))
message("Number of SNPs: ", ncol(snp))
message("Number of associations: ", ncol(cpg)*ncol(snp))
message("\n\n")


# ASSOC 
message("Running analysis")

assoc <- expand.grid(CpG = colnames(cpg), inversion = colnames(snp), stringsAsFactors = FALSE)
l <- mclapply(1:nrow(assoc), function(row){
  if(row %% 1000 == 0) message(row, " of ", nrow(assoc))
  probe <- assoc[row, "CpG"]
  inv <- assoc[row, "inversion"]
  linear_regression(cpg[, probe], snp[, inv])
}, mc.cores = threads)
l <- do.call(rbind, l)
assoc$BETA <- l[,1]
assoc$SE <- l[,2]
assoc$PVAL <- l[,3]
assoc$N <- l[,4]


message("Creating GWAMA format")

assoc <- dplyr::select(assoc, PHENOTYPE=CpG, MARKERNAME=inversion, BETA, SE, PVAL, N)
assoc <- merge(assoc, snpinfo, by="MARKERNAME", all.x=TRUE)
assoc$SNP <- assoc$MARKERNAME
assoc$MARKERNAME <- paste0(assoc$MARKERNAME, "_", assoc$PHENOTYPE)

out <- gzfile(out_file, "w")
write.table(assoc, file=out, row=FALSE, col=TRUE, qu=FALSE)
close(out)

message("Done")
