suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))

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
assoc_file <- arguments[3]
freq_file <- arguments[4]
out_file <- arguments[5]


snp <- fread(paste0("zcat ", snp_file))
cpg <- fread(cpg_file)
assoc <- fread(paste0("zcat ", assoc_file))
freq <- fread(paste0("zcat ", freq_file))

# SNP

message("Organising SNP data")

iid <- snp$IID
snp <- snp[, c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"):=NULL]

snp <- as.matrix(snp)
rownames(snp) <- iid
snpinfo <- do.call(rbind, strsplit(colnames(snp), split="_"))
snpinfo <- data.frame(SNP=snpinfo[,1], EA=snpinfo[,2], stringsAsFactors=FALSE)

colnames(snp) <- snpinfo$SNP
snpinfo <- merge(snpinfo, freq, by.x="SNP", by.y="V2", all.x=TRUE)
stopifnot(all(snpinfo$EA == snpinfo$V3))
snpinfo <- dplyr::select(snpinfo, MARKERNAME=SNP, EA=EA, NEA=V4, EAF=V5)


# CPG

message("Organising CpG data")

iid <- cpg$IID
cpg <- cpg[, c("FID", "IID"):=NULL]
cpg <- as.matrix(cpg)
rownames(cpg) <- iid

snp <- snp[rownames(snp) %in% rownames(cpg), ]
stopifnot(all(rownames(cpg) %in% rownames(snp)))
stopifnot(all(rownames(snp) %in% rownames(cpg)))

cpg <- cpg[match(rownames(snp), rownames(cpg)), ]

stopifnot(all(rownames(cpg) == rownames(snp)))


# ASSOC 

message("Organising association lists")

assoc <- subset(assoc, V2 %in% colnames(snp) & V3 %in% colnames(cpg))
assoc$x <- match(assoc$V2, colnames(snp))
assoc$y <- match(assoc$V3, colnames(cpg))

all(colnames(snp)[assoc$x] == assoc$V2)
all(colnames(cpg)[assoc$y] == assoc$V3)



# ASSOC for CIS

message("Running cis analysis")

assoc_cis <- subset(assoc, V4=="c")
l <- list()
for(i in 1:nrow(assoc_cis))
{
	if(i %% 1000 == 0) message(i, " of ", nrow(assoc_cis))
	l[[i]] <- linear_regression(cpg[, assoc_cis$y[i]], snp[, assoc_cis$x[i]])
}
l <- do.call(rbind, l)
assoc_cis$BETA <- l[,1]
assoc_cis$SE <- l[,2]
assoc_cis$PVAL <- l[,3]
assoc_cis$N <- l[,4]



# Get the best cis SNP for each trans CPG

message("Identifying best cis-SNPs for each CpG")

temp <- group_by(assoc_cis, V3) %>%
slice(which.min(PVAL)) %>%
dplyr::select(V3, cissnp=V2, cisx=x)

assoc_trans <- subset(assoc, V4=="t")
assoc_trans <- merge(assoc_trans, temp, by="V3", all.x=TRUE)

assoc_ct <- subset(assoc_trans, !is.na(cisx))
assoc_trans <- subset(assoc_trans, is.na(cisx), select=-c(cissnp, cisx))



# ASSOC for TRANS

message("Running trans analysis")

l <- list()
for(i in 1:nrow(assoc_trans))
{
	if(i %% 1000 == 0) message(i, " of ", nrow(assoc_trans))
	l[[i]] <- linear_regression(cpg[, assoc_trans$y[i]], snp[, assoc_trans$x[i]])
}
l <- do.call(rbind, l)
assoc_trans$BETA <- l[,1]
assoc_trans$SE <- l[,2]
assoc_trans$PVAL <- l[,3]
assoc_trans$N <- l[,4]



# ASSOC for TRANS adjusting for CIS

message("Running trans analysis adjusting for cis-SNPs")

l <- list()
for(i in 1:nrow(assoc_ct))
{
	if(i %% 1000 == 0) message(i, " of ", nrow(assoc_ct))
	mod <- summary(lm(cpg[, assoc_ct$y[i]] ~ snp[, assoc_ct$x[i]] + snp[, assoc_ct$cisx[i]]))
	l[[i]] <- c(coefficients(mod)[2,c(1,2,4)], mod$df[2])
}
l <- do.call(rbind, l)
assoc_ct$BETA <- l[,1]
assoc_ct$SE <- l[,2]
assoc_ct$PVAL <- l[,3]
assoc_ct$N <- l[,4]


message("Creating GWAMA format")

assoc <- rbind(assoc_ct, assoc_cis, assoc_trans, fill=TRUE)
assoc <- dplyr::select(assoc, PHENOTYPE=V3, MARKERNAME=V2, CISTRANS=V4, BETA, SE, PVAL, N)
assoc <- merge(assoc, snpinfo, by="MARKERNAME", all.x=TRUE)
assoc$SNP <- assoc$MARKERNAME
assoc$MARKERNAME <- paste0(assoc$MARKERNAME, "_", assoc$PHENOTYPE)

out <- gzfile(out_file, "w")
write.table(assoc, file=out, row=FALSE, col=TRUE, qu=FALSE)
close(out)

message("Done")
