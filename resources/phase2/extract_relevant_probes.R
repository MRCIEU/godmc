arguments <- commandArgs(T)

methylationfile <- arguments[1]
fam_file <- arguments[2]
probedir <- arguments[3]
phase2_beta_prefix <- arguments[4]


message("Reading methylation data...")
load(methylationfile)


message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")


fam <- read.table(fam_file, he=F, stringsAsFactors = FALSE)
names(fam)[1:2] <- c("FID","IID")
rownames(fam) <- fam$IID
fam <- subset(fam, IID %in% colnames(norm.beta), select=-c(V3,V4,V5,V6))
norm.beta <- norm.beta[, colnames(norm.beta) %in% fam$IID]
index <- match(colnames(norm.beta), fam$IID)
fam <- fam[index,]

# Make sure that fam IID is the same order as norm.beta IDs
stopifnot(all(fam$IID == colnames(norm.beta)))
fam <- as.matrix(fam)


# Get the probe lists

file_list <- list.files(probedir, pattern=".probes")
n_sets <- length(file_list)

message("Extracting CpGs from ", n_sets," probe sets")

for (i in 1:n_sets)
{
	message(i, " of ", n_sets)
	probes <- scan(file.path(probedir, file_list[i]), what="character")
	m <- t(norm.beta[rownames(norm.beta) %in% probes, ])
	nom <- colnames(m)
	m <- matrix(as.character(m), nrow(m), ncol(m))
	colnames(m) <- nom
	m <- cbind(fam, m)

	out <- paste0(phase2_beta_prefix, i)

	write.table(m, file=out, row=FALSE, col=TRUE, quote=FALSE)
}

message("Successfully extracted probes from beta matrix")
 
