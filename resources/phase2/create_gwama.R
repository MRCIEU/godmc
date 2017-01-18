library(data.table)

arguments <- commandArgs(T)
results_file <- arguments[1]
freq_file <- arguments[2]

results <- fread(results_file)
results$index <- 1:nrow(results)
names(results) <- c("PHENOTYPE", "MARKERNAME", "N", "BETA", "SE", "PVAL", "CISTRANS", "index")


freq <- fread(paste0("zcat ", freq_file))
freq <- subset(freq, select=c(V2, V3, V4, V5))
names(freq) <- c("MARKERNAME", "EA", "NEA", "EAF")

dat <- merge(results, freq, by="MARKERNAME", all.x=TRUE)
dat <- dat[order(dat$index), ]
dat <- subset(dat, select=c(PHENOTYPE, CISTRANS, MARKERNAME, EA, NEA, BETA, SE, N, EAF, PVAL))

write.table(dat, file=results_file, row=FALSE, col=TRUE, qu=FALSE)
