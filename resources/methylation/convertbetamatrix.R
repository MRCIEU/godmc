library(parallel)
suppressMessages(library(matrixStats))

	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	no.probesets <- as.numeric(arguments[2])
    methylationdir<-arguments[3]
    probedir <- arguments[4]
    fam_file <- arguments[5]

    message("Reading methylation data...")
	load(methylationfile)
    message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")
    
    fam <- read.table(fam_file, he=F)
	names(fam)[1:2]<-c("FID","IID")
	rownames(fam)<-fam$IID
	fam <- subset(fam, IID %in% colnames(norm.beta), select=-c(V3,V4,V5,V6))
    norm.beta <- norm.beta[, colnames(norm.beta) %in% fam$IID]
    
   
    message("Extracted " ,no.probesets," probe subsets")
    
    probefiles=list.files(path=probedir,pattern=".allcohorts.probes")
    
    message("Found " ,length(probefiles)," probefiles")

    for (p in 1:no.probesets){
    
    probes <- read.table(paste(probedir,probefiles[p],sep="/"), he=F, stringsAsFactors=F)
    p1<-gsub("cis_trans.","",probefiles[p])
    p1<-gsub(".allcohorts.probes","",p1)
    index <- apply(probes, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
    probes <- probes[!index, ]
    probes<-probes[probes%in%row.names(norm.beta)]
    m<-match(probes,rownames(norm.beta))
    norm.beta.subset <- norm.beta[m,]
    norm.beta.subset <- t(norm.beta.subset)
    
    m<-match(fam$IID,rownames(norm.beta.subset))
    norm.beta.subset<-data.frame(fam,norm.beta.subset[m,])

    #index <- which(is.na(norm.beta), arr.ind = TRUE)
	#if (length(index)>0){
    #message("Replace ",length(index)," missing values with rowmeans")
    #norm.beta[index] <- rowMeans(norm.beta, na.rm = TRUE)[index[, "row"]] }

	write.table(norm.beta.subset, paste(methylationdir,"/plink.methylation.subset.",p1,".txt",sep=""),sep=" ",na = "NA",col.names=T,row.names=F,quote=F)
	message("Successfully extracted probes from beta matrix for probeset ",p)
    }
    
     

    
    


