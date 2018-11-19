args <- commandArgs(trailingOnly=TRUE)
pheno_file <- args[1]
pheno_file2 <- args[2]
pheno <- args[3]
out <- args[4]

phenodf <- read.table(pheno_file, header = TRUE)
phenoids <- read.table(pheno_file2, header = FALSE)

resdf <- data.frame(FID = phenoids[, 1], IID = phenoids[, 2], phenodf[, pheno])
write.table(resdf, col.names = FALSE, row.names = FALSE, file = paste0(out, pheno, ".plink"),
            quote = FALSE)