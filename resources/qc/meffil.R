arguments <- commandArgs(T)
idatdir <- arguments[1]
geno_probesnps <- arguments[2]
naeemfile <- arguments[3]
famfile <- arguments[3]
savefile <- arguments[4]


library(meffil)
options(mc.cores=16)

# Get control genotypes
genotypes <- meffil.extract.genotypes(geno_probesnps)

# Get sex from famfile
fam <- read.table(famfile, he=F, stringsAsFactors=FALSE)
sex <- subset(fam, select=c(V2, V5))
names(sex) <- c("IID", "sex")
sex$sex[sex==0] <- "F"
sex$sex[sex==1] <- "M"

# Generate samplesheet
# This needs to be done manually
# The function below will try to guess the sample information from the idat filenames
# Make sure that the Sample_Name column corresponds to the IID column in genotypes
# Make sure that the Sex column is completed with "M" and "F" entries.
samplesheet <- meffil.create.samplesheet("~/data/test_meffil/")


# Background and dye bias correction, sexprediction, cell counts estimates
qc.objects <- meffil.qc(samplesheet, cell.type.reference="blood gse35069", verbose=TRUE)


# Generate QC report
qc.summary <- meffil.qc.summary(qc.objects, genotypes=genotypes)
meffil.qc.report(qc.summary, output.file="qc-report.html")

# Remove outlier samples if necessary
qc.objects <- meffil.remove.samples(qc.objects, qc.summary$bad.samples$sample.name)

# Perform quantile normalization
norm.objects <- meffil.normalize.quantiles(qc.objects, random.effects="Slide", number.pcs=10)

# Generate normalized probe values
naeemfilter <- subset(naeem, Flag.discard.keep. == "discard")$probe
cpgs_to_remove <- unique(c(naeemfilter, qc.summary$bad.cpgs$name)
norm.beta <- meffil.normalize.samples(norm.objects, cpglist.remove=cpgs_to_remove)

# Generate normalization report
norm.summary <- meffil.normalization.summary(norm.beta, norm.objects)
meffil.normalization.report(norm.summary, output.file="normalization-report.html")