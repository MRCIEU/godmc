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
# FID IID Sex Smoking

if("Sex" %in% names(allcovs))
{
	allcovs$Sex[allcovs$Sex == "M"] <- 1
	allcovs$Sex[allcovs$Sex == "F"] <- 2
	allcovs$Sex <- as.numeric(allcovs$Sex)
}


if("Sex" %in% names(allcovs))
{
	dat <- merge(fam, subset(allcovs, select=c("IID", "Sex")), by.x="V2", by.y="IID")
	dat <- merge(dat, smok, by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", "Sex", "Smoking"))
	write.table(dat, paste0(out_file, ".aar"), row=F, col=F, qu=F)
}

# Create covariates for Smoking GWAS
# FID IID Age Sex

covnames <- c("Age", "Sex")
covnames <- covnames[covnames %in% names(allcovs)]
if(length(covnames) > 0)
{
	dat <- merge(fam, subset(allcovs, select=c("IID", covnames)), by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", covnames))
	write.table(dat, paste0(out_file, ".smoking"), row=F, col=F, qu=F)
}


# Create covariates for cellcounts GWAS
# FID IID Smoking Age Sex

covnames <- c("Age", "Sex")
covnames <- covnames[covnames %in% names(allcovs)]
dat <- merge(fam, smok, by.x="V2", by.y="IID")
dat <- subset(dat, select=c("V1", "V2", "Smoking"))
if(length(covnames) > 0)
{
	dat <- merge(dat, subset(allcovs, select=c("IID", covnames)), by.x="V2", by.y="IID")
	dat <- subset(dat, select=c("V1", "V2", "Smoking", covnames))
}
write.table(dat, paste0(out_file, ".cellcounts"), row=F, col=F, qu=F)

