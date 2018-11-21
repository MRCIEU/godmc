#'#################################################################################
#'#################################################################################
#' Run scoreInvHap in GoDMC samples
#'#################################################################################
#'#################################################################################
args <- commandArgs(trailingOnly=TRUE)
plink <- args[1]
out_folder <- args[2]
res_folder <- args[3]
threads <- args[4]

library(scoreInvHap)
library(snpStats)
library(BiocParallel)

data(inversionGR)
invs <- names(inversionGR)

## Remove inversion in Chr X
invs <- invs[invs != "invX_006"]

snps <- read.plink(plink, na.strings = "")
snps <- checkSNPs(snps)$genos

invstat<-lapply(invs, function (inv) {
  print(inv)
  result <- tryCatch(scoreInvHap(snps, BPPARAM = MulticoreParam(threads), inv=inv),
                     error = function(e) as.character(e))
                     }
)
names(invstat) <- invs

## Extract inversions with genotyping errors
errors <- unlist(invstat[sapply(invstat, class) == "character"])
errorInvs <- names(errors)

invstat <- invstat[sapply(invstat, class) == "scoreInvHapRes"]
invstat2 <- invstat[sapply(invstat, function(x) max(numSNPs(x))) > 3 ]
snpInvs <- names(invstat)[!names(invstat) %in% names(invstat2)]
  
#Crear data.frame
newClassification <- function(object, minDiff = 0){
  scores <- scores(object)
  scores[scores == Inf] <- 1
  
  scores2 <- data.frame(NN = apply(scores[, grep("N.*N", colnames(scores)), drop = FALSE], 1, max),
                        NI = apply(scores[, grep("N.*I", colnames(scores)), drop = FALSE], 1, max), 
                        II = apply(scores[, grep("I.*I", colnames(scores)), drop = FALSE], 1, max))

  class <- colnames(scores2)[max.col(scores2)]
  names(class) <- rownames(scores2)
  
  diff <- apply(scores2, 1, function(x) max(x) - sort(x, decreasing = TRUE)[2])
  return(class[diff > minDiff])
}
invsList <- lapply(invstat2, function(x) newClassification(x))

## Remove inversions with more 10% without genotype
invsList <- invsList[lengths(invsList) > 0.9*max(lengths(invsList))]
ids <- rownames(snps$genotypes)

genoInvs <- names(invstat2)[!names(invstat2) %in% names(invsList)]

invsDF <- Reduce(function(x, y) rbind(x, y[ids]), invsList)
invsDF <- cbind(id = names(invsList), invsDF)
invsDF[invsDF == "NN"] <- "0"
invsDF[invsDF == "NI"] <- "1"
invsDF[invsDF == "II"] <- "2"
write.table(invsDF, file = paste0(out_folder, "/inversionQTL.txt"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

rownames(invsDF) <- invsDF[, 1]
invsDF <- invsDF[, -1]

class(invsDF) <- "numeric"
mat <- new("SnpMatrix", t(invsDF)+1)
annot <- as.data.frame(inversionGR[colnames(mat)])
annot$genetic.distance <- 0
annot$allele.1 <- "N"
annot$allele.2 <- "I"
annot$chr <- gsub("chr", "", as.character(annot$seqnames))

write.plink(snps = mat, file.base = paste0(out_folder, "/inversions"), na.code = NA,
            subject.data = snps$fam,
            pedigree =   as.character(pedigree),
            id = as.character(member), father = as.character(father), 
            mother = as.character(mother), sex = as.numeric(sex), 
            phenotype = as.numeric(affected),
            chromosome = annot$chr, position = annot$start,
            allele.1 = annot$allele.1, allele.2 = annot$allele.2, 
            genetic.distance = annot$genetic.distance)

df <- data.frame(inversions = c(errorInvs, snpInvs, genoInvs), 
                 Reason = c(errors, rep(c("Not enough SNPs", "Bad genotyping"), 
                              c(length(snpInvs), length(genoInvs)))))
write.table(df, file = paste0(res_folder, "/badInversions.txt"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)

invFreqs <- round(rowMeans(invsDF, na.rm = TRUE)/2*100, 2)
df <- data.frame(inversion = names(invFreqs), freq = invFreqs)
write.table(df, file = paste0(res_folder, "/inversionsSummary.txt"), sep = "\t",
            row.names = FALSE, col.names = TRUE, quote = FALSE)
