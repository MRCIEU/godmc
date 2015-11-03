iids <- read.table("data.fam")
iids <- as.factor(iids[,2])

KINGmat <- king2mat(
	file.kin0 = "king.kin0", 
	file.kin = "king.kin", 
	iids = iids
)




library(SNPRelate)

bed.fn <- "data.bed"
bim.fn <- "data.bim"
fam.fn <- "data.fam"
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "data.gds")

genofile <- snpgdsOpen("data.gds")


ibd <- snpgdsIBDKING(genofile)


geno <- GdsGenotypeReader(filename="data.gds")
genoData <- GenotypeData(geno)

ibd <- ibd$kinship
colnames(ibd) <- iids
rownames(ibd) <- iids

mypcair <- pcair(genoData = genoData, kinMat = ibd, divMat = ibd)