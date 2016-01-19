write_covs <- function(dat, filename)
{
	fcols <- grep("_factor", names(dat))
	ncols <- grep("_numeric", names(dat))
	fdat <- data.frame(dat[,1:2], dat[,fcols])
	ndat <- data.frame(dat[,1:2], dat[,ncols])
	write.table(fdat, file=paste0(filename, ".factor"), row=F, col=F, qu=F)
	write.table(ndat, file=paste0(filename, ".numeric"), row=F, col=F, qu=F)
	write.table(dat, file=filename, row=F, col=F, qu=F)
}


arguments <- commandArgs(T)

covs_file <- arguments[1]
aar_file <- arguments[2]
smok_file <- arguments[3]
fam_file <- arguments[4]
out_file <- arguments[5]

allcovs <- read.table(covs_file, he=T, stringsAsFactors=FALSE)
aar <- read.table(aar_file, he=T, stringsAsFactors=FALSE)
smok <- read.table(smok_file, he=T, stringsAsFactors=FALSE)
fam <- read.table(fam_file, stringsAsFactors=FALSE)[,1:2]

# Create covariates for AAR GWAS
# FID IID Sex_factor Smoking

if("Sex_factor" %in% names(allcovs))
{
	allcovs$Sex_factor[allcovs$Sex_factor == "M"] <- 1
	allcovs$Sex_factor[allcovs$Sex_factor == "F"] <- 2
	allcovs$Sex_factor <- as.numeric(allcovs$Sex_factor)
}


if("Sex_factor" %in% names(allcovs))
{
	dat <- merge(fam, subset(allcovs, select=c("IID", "Sex_factor")), by.x="V2", by.y="IID")
	dat <- merge(dat, smok, by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", "Sex_factor", "Smoking"))
	names(dat)[names(dat) == "Smoking"] <- "Smoking_numeric"
	write_covs(dat, paste0(out_file, ".aar"))
}

# Create covariates for Smoking GWAS
# FID IID Age_numeric Sex

covnames <- c("Age_numeric", "Sex_factor")
covnames <- covnames[covnames %in% names(allcovs)]
if(length(covnames) > 0)
{
	dat <- merge(fam, subset(allcovs, select=c("IID", covnames)), by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", covnames))
} else {
	dat <- fam[,1:2]
}
write_covs(dat, paste0(out_file, ".smoking"))


# Create covariates for cellcounts GWAS
# FID IID Smoking Age_numeric Sex

covnames <- c("Age_numeric", "Sex_factor")
covnames <- covnames[covnames %in% names(allcovs)]
dat <- merge(fam, smok, by.x="V2", by.y="IID")
dat <- subset(dat, select=c("V1", "V2", "Smoking"))
if(length(covnames) > 0)
{
	dat <- merge(dat, subset(allcovs, select=c("IID", covnames)), by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", "Smoking", covnames))
}
names(dat)[names(dat) == "Smoking"] <- "Smoking_numeric"
write_covs(dat, paste0(out_file, ".cellcounts"))
