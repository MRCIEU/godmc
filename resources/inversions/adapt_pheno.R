args <- commandArgs(trailingOnly=TRUE)
pheno_file <- args[1]
pheno <- args[2]
out <- args[3]

phenodf <- read.table(pheno_file, header = TRUE)

resdf <- data.frame(A = phenodf$IID, B = phenodf$IID, Height = phenodf[, "Height"])
write.table(resdf, col.names = FALSE, row.names = FALSE, file = paste0(out, pheno, ".plink"),
            quote = FALSE)