library(parallel)
suppressMessages(library(matrixStats))

	arguments <- commandArgs(T)

	methylationfile <- arguments[1]
	cov_file <- arguments[2]
	no.probesets <- as.numeric(arguments[3])
    probedir<-arguments[4]
    methylationdir <-arguments[5]
    fam_file <- arguments[6]
    cov_out <-arguments[7]
    nthreads<-as.numeric(arguments[8])
  

    message("Reading methylation data...")
	load(methylationfile)
    message("Data size: ", ncol(norm.beta), " individuals and ", nrow(norm.beta), " CpGs.")
    

    covs <- read.table(cov_file, he=T)
	index <- apply(covs, 1, function(x) any(is.na(x) | is.nan(x) | is.infinite(x)))
	covs <- covs[!index, ]
	rownames(covs) <- covs$IID
	covs <- subset(covs, IID %in% colnames(norm.beta), select=-c(IID))
    norm.beta <- norm.beta[, colnames(norm.beta) %in% rownames(covs)]

    fam <- read.table(fam_file, he=F)
	names(fam)[1:2]<-c("FID","IID")
	rownames(fam)<-fam$IID
	fam <- subset(fam, IID %in% colnames(norm.beta), select=-c(V3,V4,V5,V6))
    norm.beta <- norm.beta[, colnames(norm.beta) %in% fam$IID]
    
    message("Identifying methylation outliers")

	niter <- 3
	outlier_threshold <- 10
	norm.beta.copy <- norm.beta
	for(i in 1:niter)
	{
		sds <- rowSds(norm.beta.copy, na.rm=T)
		means <- rowMeans(norm.beta.copy, na.rm=T)
		norm.beta.copy[norm.beta.copy > means + sds*outlier_threshold | norm.beta.copy < means - sds*outlier_threshold] <- NA
	}
	outlier_count <- apply(norm.beta.copy, 1, function(x) sum(is.na(x)))
	norm.beta.copy <- is.na(norm.beta.copy)
    norm.beta[norm.beta.copy] <- NA

    message("Rank transform methylation betas")

    norm.beta2<-mclapply(1:dim(norm.beta)[1], rntransform,mc.cores=nthreads)
    norm.beta<-norm.beta2

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

    write.table(norm.beta.subset, paste(methylationdir,"/gcta.methylation.subset.",p1,".txt",sep=""),sep=" ",na = "NA",col.names=T,row.names=F,quote=F)
    message("Successfully extracted probes from beta matrix for probeset ",p)
    }
    
    g1<-grep("genetic_pc",names(covs))
    covs<-covs[,-g1]
    
    g<-grep("factor",names(covs))
    
    covs.f<-covs[,g]
    m<-match(fam$IID,rownames(covs.f))
    covs.f<-data.frame(fam,covs.f[m,])
    write.table(covs.f,paste(cov_out,".factor",sep="") ,sep=" ",na = "NA",col.names=F,row.names=F,quote=F)    
   
    covs.n<-covs[,-g]
    m<-match(fam$IID,rownames(covs.n))
    covs.n<-data.frame(fam,covs.n[m,])
    write.table(covs.n,paste(cov_out,".numeric",sep="") ,sep=" ",na = "NA",col.names=F,row.names=F,quote=F)    
   
    rntransform <- function(x)
{
    out <- rank(x) - 0.5
    out[is.na(x)] <- NA
    mP <- 0.5/max(out, na.rm = T)
    out <- out/(max(out, na.rm = T) + 0.5)
    out <- scale(qnorm(out))
    out
}

    
    


