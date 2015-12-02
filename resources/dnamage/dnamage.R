## USE:

## http://labs.genetics.ucla.edu/horvath/dnamage/TUTORIAL1.pdf
## http://labs.genetics.ucla.edu/horvath/dnamage/probeAnnotation21kdatMethUsed.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/datMiniAnnotation27k.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/AdditionalFile3.csv
## http://labs.genetics.ucla.edu/horvath/dnamage/StepwiseAnalysis.txt
## http://labs.genetics.ucla.edu/horvath/dnamage/NORMALIZATION.R

# library(RPMM)
# library(sqldf) 
## http://cran.r-project.org/src/contrib/Archive/sqldf/sqldf_0.4-6.4.tar.gz
## http://cran.r-project.org/src/contrib/Archive/RSQLite.extfuns/RSQLite.extfuns_0.0.1.tar.gz
# library(impute) ## from bioconductor
# library(WGCNA)


### Normalisation functions
### ORIGINAL AUTHOR: Andrew Teschendorff
# The original BMIQ function from Teschendorff 2013 adjusts for the type-2 bias in 
# Illumina Infinium 450k data.
# Later functions and edits were provided by yours truly, Steve Horvath.
# I changed the code so that one can calibrate methylation data to a gold standard.
# Specifically, I took version v_1.2 by Teschendorff  and fixed minor issues. 
# Also I made the code more robust e.g. by changing the optimization algorithm.
# Toward this end, I used the method="Nelder-Mead" in optim()

### Later functions and edits by Steve Horvath 
### # Steve Horvath took version v_1.2 by Teschendorff 
# and fixed minor errors. Also he made the code more robust.
# Importantly, SH changed the optimization algorithm to make it #more robust. 
# SH used method="Nelder-Mead" in optim() since the other #optimization method sometimes gets stuck.
#Toward this end, the function blc was replaced by blc2.




main <- function()
{
	# Args
	arguments <- commandArgs(T)
	beta_file <- arguments[1]
	cov_file <- arguments[2]
	fam_file <- arguments[3]
	out_file <- arguments[4]
    age.acceleration.plot<-arguments[5]
    SD<-as.numeric(arguments[6])

	message("Getting probe parameters")
	dnamage.probeAnnotation21kdatMethUsed=read.csv("resources/dnamage/probeAnnotation21kdatMethUsed.csv.gz")
	dnamage.probeAnnotation27k=read.csv("resources/dnamage/datMiniAnnotation27k.csv.gz")
	dnamage.datClock=read.csv("resources/dnamage/AdditionalFile3.csv.gz")

	message("Reading in data")
	covs <- read.table(cov_file, header=T, stringsAsFactors=FALSE)
	fam <- read.table(fam_file, stringsAsFactors=FALSE)
	m<-match(fam[,2],covs[,1])
	covs<-covs[m,]

	index <- grep("^age_numeric$", names(covs), ignore.case=TRUE)
	if(length(index) != 1)
	{
		stop("There should be only one column in the covariate file called 'Age_numeric', regardless of case.")
	}
	names(covs)[index] <- "Age_numeric"
	if("FID" %in% names(covs)) covs <- subset(covs, select=-c(FID))
	load(beta_file)
    m<-match(fam[,2],colnames(norm.beta))
    norm.beta<-norm.beta[,m]
	message("Predicting age")
	agepred <- dnamage(norm.beta, normalizeData=FALSE, dnamage.probeAnnotation27k, dnamage.probeAnnotation21kdatMethUsed, dnamage.datClock)
	agepred$SampleID <- colnames(norm.beta)


	message("Generating age accelerated residuals")
	phen <- generate.aar(agepred, covs)

	pdf(age.acceleration.plot, width=12, height=8)
	par(mfrow=c(2,2))
	plot(phen$Age_numeric,phen$DNAmAge,cex.main=0.7,cex=0.7, main=paste("correlation between predicted age and actual age=",cor(phen$Age_numeric, phen$DNAmAge,use="pair"),sep=""))

    write.table(phen, file=paste0(out_file, ".txt"), row=F, col=T, qu=F)

	nom <- names(phen)
	#phen <- merge(phen, fam, by.x="IID", by.y="V2")
	#phen <- subset(phen, select=c(V1, IID, AAR))
	m<- match(fam$V2,phen$IID)
	phen<-data.frame(FID=fam$V1,IID=phen$IID[m],AAR=phen$AAR[m])
	write.table(phen, file=paste0(out_file, ".aar.plink"), row=F, col=F, qu=F)


	par(mfrow=c(2,2))
	plot(phen$AAR, xlab="", main=paste("age acceleration (N=", length(which(!is.na(phen$AAR))),")",sep=""),cex.main=0.7)
	hist(phen$AAR, xlab="", main=paste("age acceleration (N=", length(which(!is.na(phen$AAR))),")",sep=""),cex.main=0.7)
	abline(v=mean(phen$AAR,na.rm=T)-SD*sd(phen$AAR,na.rm=T),lty=2)
	abline(v=mean(phen$AAR,na.rm=T)+SD*sd(phen$AAR,na.rm=T),lty=2)
	qqnorm(phen$AAR, main=paste("age acceleration (N=", length(which(!is.na(phen$AAR))),"; shapiroP=",signif(as.numeric(shapiro.test(phen$AAR)[2]),2),")",sep=""),cex.main=0.7)
	qqline(phen$AAR)
	quiet <- dev.off()
}



dnamage <- function(x,normalizeData=TRUE, dnamage.probeAnnotation27k, dnamage.probeAnnotation21kdatMethUsed, dnamage.datClock) {
	dat0 <- data.frame(id=rownames(x), x)
	nSamples=dim(dat0)[[2]]-1
	nProbes= dim(dat0)[[1]]
	## the following command may not be needed. But it is sometimes useful when you use read.csv.sql
	dat0[,1]= gsub(x=dat0 [,1],pattern="\"",replacement="")
	cat(paste( "The methylation data set contains", nSamples, "samples (e.g. arrays) and ",
			  nProbes, " probes.\n"))
	if (nSamples==0) {
		stop("ERROR: There must be a data input error since there seem to be no samples.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql(). Samples correspond to columns in that file .")
	}
	if (nProbes==0) {
		stop("ERROR: There must be a data input error since there seem to be zero probes.\n Make sure that you input a comma delimited file (.csv file)\n that can be read using the R command read.csv.sql CpGs correspond to rows.")
	}
	if ( nSamples > nProbes ) {
		warning("MAJOR WARNING: It worries me a lot that there are more samples than CpG probes.\n Make sure that probes correspond to rows and samples to columns.\n I wonder whether you want to first transpose the data and then resubmit them? In any event, I will proceed with the analysis.")
	}
	if ( is.numeric(dat0[,1]) ) {
		stop(paste("Error: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]))
	}
	if ( !is.character(dat0[,1]) ) {
		warning(paste("Major Warning: The first column does not seem to contain probe identifiers (cg numbers from Illumina) since these entries are numeric values. Make sure that the first column of the file contains CpG probe identifiers such as cg00000292. Instead it contains ", dat0[1:3,1]))
	}


	nonNumericColumn=rep(FALSE, dim(dat0)[[2]]-1)
	for (i in 2:dim(dat0)[[2]] ){ nonNumericColumn[i-1]=! is.numeric(dat0[,i]) }
	if ( sum(nonNumericColumn) >0 ) {
		warning(paste("MAJOR WARNING: Possible input error. The following samples contain nonnumeric beta values: ", colnames(dat0)[-1][ nonNumericColumn], "\n Hint: Maybe you use the wrong symbols for missing data. Make sure to code missing values as NA in the Excel file. To proceed, I will force the entries into numeric values but make sure this makes sense."))
	}



	XchromosomalCpGs=as.character(dnamage.probeAnnotation27k$Name[dnamage.probeAnnotation27k$Chr=="X"])
	selectXchromosome=is.element(dat0[,1], XchromosomalCpGs )
	selectXchromosome[is.na(selectXchromosome)]=FALSE
	meanXchromosome=rep(NA, dim(dat0)[[2]]-1)
	if ( sum(selectXchromosome) >=500 ) {
		meanXchromosome= as.numeric(apply( as.matrix(dat0[selectXchromosome,-1]),2,mean,na.rm=TRUE))
	}

	if ( sum(is.na(meanXchromosome)) >0 ) {
		warning("Comment: There are lots of missing values for X chromosomal probes for some of the samples. This is not a problem when it comes to estimating age but I cannot predict the gender of these samples.")
		match1=match(dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])
		if ( sum( is.na(match1))>0 ) {
			missingProbes= dnamage.probeAnnotation21kdatMethUsed$Name[!is.element( dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])]
			##MATT used to stop on this error, changed to a warning
			warning(paste("Input error: You forgot to include the following ", length(missingProbes), " CpG probes (or probe names):\n "))#MATT##, paste( missingProbes, sep="",collapse=", ")))
		}
	}
	
	match1=match(dnamage.probeAnnotation21kdatMethUsed$Name , dat0[,1])
	if ( sum( is.na(match1))>0 ) ##MATT used to stop on this error, changed to a warning
		warning(paste(sum( is.na(match1)), "CpG probes cannot be matched"))
	dat1= dat0[match1,]
	asnumeric1=function(x) {as.numeric(as.character(x))}
	dat1[,-1]=apply(as.matrix(dat1[,-1]),2,asnumeric1)
	
	set.seed(1)
	## Do you want to normalize the data (recommended)?
	##MATT##normalizeData=TRUE

	##MATT added this so that all age CpGs are at least included (check for all CpGs above is a warning now)
	selectCpGsClock=is.element(dat0[,1], as.character(dnamage.datClock$CpGmarker[-1]))
	if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {stop("The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.")}
	if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {stop("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).")}


	##MATT##source("dnamage/StepwiseAnalysis.txt")
	dnamage.stepwise(dat1, meanXchromosome,
					 normalizeData=normalizeData, dnamage.probeAnnotation27k, dnamage.probeAnnotation21kdatMethUsed, dnamage.datClock)
}



dnamage.stepwise <- function(dat1, meanXchromosome,
							 normalizeData, dnamage.probeAnnotation27k, dnamage.probeAnnotation21kdatMethUsed, dnamage.datClock) {

	nSamples=dim(dat1)[[2]]-1
	nProbes= dim(dat1)[[1]]


	trafo <- function(x,adult.age=20) { 
		x=(x+1)/(1+adult.age);
		y=ifelse(x<=1, log( x),x-1);
		y
	}
	anti.trafo <- function(x,adult.age=20) {
		ifelse(x<0,
			   (1+adult.age)*exp(x)-1,
			   (1+adult.age)*x+adult.age)
	}

	## Steve Horvath: Estimating DNAm age.
	## This file assumes a data frame exists called dat1 whose rows correspond to CpGs
	## and whose first column reports the CpG identifier
	## and whose remaining columns corresponds to samples (e.g. Illumina arrays).
	
	fastImputation=FALSE
	
	##STEP 1: DEFINE QUALITY METRICS
	
	meanMethBySample =as.numeric(apply(as.matrix(dat1[,-1]),2,mean,na.rm=TRUE))
	minMethBySample   =as.numeric(apply(as.matrix(dat1[,-1]),2,min,na.rm=TRUE))
	maxMethBySample  =as.numeric(apply(as.matrix(dat1[,-1]),2,max,na.rm=TRUE))
	
	datMethUsed= t(dat1[,-1])
	colnames(datMethUsed)=as.character(dat1[,1])
	
	noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
	# print(table(noMissingPerSample))
	
	##STEP 2: Imputing 
	if (! fastImputation & nSamples>1 & max(noMissingPerSample,na.rm=TRUE)<3000 ){
		
		## run the following code if there is at least one missing
		if ( max(noMissingPerSample,na.rm=TRUE)>0 ){
			dimnames1=dimnames(datMethUsed)
			datMethUsed= data.frame(t(impute.knn(t(datMethUsed))$data))
			dimnames(datMethUsed)=dimnames1
		} # end of if
	} # end of if (! fastImputation )
	
	if ( max(noMissingPerSample,na.rm=TRUE)>=3000 ) fastImputation=TRUE
	
	if ( fastImputation | nSamples==1 ){
		noMissingPerSample=apply(as.matrix(is.na(datMethUsed)),1,sum)
		table(noMissingPerSample)
		if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) >= 3000 ) {normalizeData=FALSE}

		## run the following code if there is at least one missing
		if ( max(noMissingPerSample,na.rm=TRUE)>0 & max(noMissingPerSample,na.rm=TRUE) < 3000 ){
			dimnames1=dimnames(datMethUsed)
			for (i in which(noMissingPerSample>0) ){
				selectMissing1=is.na(datMethUsed[i,])
				datMethUsed[i,selectMissing1] = as.numeric(dnamage.probeAnnotation21kdatMethUsed$goldstandard2[selectMissing1])
			} # end of for loop
			dimnames(datMethUsed)=dimnames1
		} # end of if
	} # end of if (! fastImputation )


	## STEP 3: Data normalization (each sample requires about 8 seconds). It would be straightforward to parallelize this operation.
	
	if (normalizeData ){
		datMethUsedNormalized=BMIQcalibration(datM=datMethUsed,goldstandard.beta= dnamage.probeAnnotation21kdatMethUsed$goldstandard2,plots=FALSE)
	}
	if (!normalizeData ){ datMethUsedNormalized=datMethUsed }
	rm(datMethUsed); gc()
	
	##STEP 4: Predict age and create a data frame for the output (referred to as datout)
	selectCpGsClock=is.element(dimnames(datMethUsedNormalized)[[2]], as.character(dnamage.datClock$CpGmarker[-1]))
	if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {stop("The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.")}
	if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {stop("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).")}
	if (nSamples>1 ) {
		datMethClock0=data.frame(datMethUsedNormalized[,selectCpGsClock])
		datMethClock= data.frame(datMethClock0[ as.character(dnamage.datClock$CpGmarker[-1])])
		dim(datMethClock)
		predictedAge=as.numeric(anti.trafo(dnamage.datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dnamage.datClock$CoefficientTraining[-1])))
	} # end of if
	
	
	if (nSamples==1 ) {
		datMethUsedNormalized2=data.frame(rbind(datMethUsedNormalized,datMethUsedNormalized))
		datMethClock0=data.frame(datMethUsedNormalized2[,selectCpGsClock])
		datMethClock= data.frame(datMethClock0[ as.character(dnamage.datClock$CpGmarker[-1])])
		dim(datMethClock)
		predictedAge=as.numeric(anti.trafo(dnamage.datClock$CoefficientTraining[1]+as.matrix(datMethClock)%*% as.numeric(dnamage.datClock$CoefficientTraining[-1])))
		predictedAge=predictedAge[1]
	} # end of if
	
	
	
	## Let's add comments to the age prediction
	Comment=ifelse ( predictedAge <0, "Negative DNAm age.", ifelse ( predictedAge >100, "Old DNAm age.", rep("",length(predictedAge))))
	
	Comment[is.na(predictedAge)]="Age prediction was not possible. "
	
	
	if ( sum( selectCpGsClock) < dim(dnamage.datClock)[[1]]-1 ) {
		Comment=rep("ERROR: The CpGs listed in column 1 of the input data did not contain the CpGs needed for calculating DNAm age. Make sure to input cg numbers such as cg00075967.",length(predictedAge) )}
	
	
	if ( sum( selectCpGsClock) > dim(dnamage.datClock)[[1]]-1 ) {
		Comment=rep("ERROR: The CpGs listed in column 1 of the input data contain duplicate CpGs. Each row should report only one unique CpG marker (cg number).",length(predictedAge) )}
	
	
	restSamples=-minMethBySample>0.05 | maxMethBySample>1.05;
	restSamples[is.na(restSamples)]=FALSE
	lab1="MAJOR WARNING: Probably you did not input beta values since either minMethBySample<-0.05 or maxMethBySample>1.05.";Comment[restSamples]= paste(Comment[restSamples],lab1)
	
	restSamples= noMissingPerSample >0 & noMissingPerSample <=100;lab1="WARNING: Some beta values were missing, see noMissingPerSample."; Comment[restSamples]= paste(Comment[restSamples],lab1)
	restSamples= noMissingPerSample >3000;lab1="MAJOR WARNING: More than 3k missing values!!"; Comment[restSamples]= paste(Comment[restSamples],lab1)
	
	restSamples= noMissingPerSample >100 & noMissingPerSample <=3000 ;lab1="MAJOR WARNING: noMissingPerSample>100"
	Comment[restSamples]= paste(Comment[restSamples],lab1)
	restSamples=meanMethBySample>.35;
	restSamples[is.na(restSamples)]=FALSE
	lab1="Warning: meanMethBySample is >0.35";Comment[restSamples]= paste(Comment[restSamples],lab1)
	restSamples=meanMethBySample<.25;
	restSamples[is.na(restSamples)]=FALSE; lab1="Warning: meanMethBySample is <0.25"
	Comment[restSamples]= paste(Comment[restSamples],lab1)
	datout=data.frame(SampleID=colnames(dat1)[-1], DNAmAge=predictedAge, Comment, noMissingPerSample,meanMethBySample, minMethBySample, maxMethBySample)
	
	# The thresholds to predict gender are incorrect.
	#if ( !is.null( meanXchromosome) ){  
		
	#	if ( length( meanXchromosome)==dim(datout)[[1]] ){
#			predictedGender=ifelse(meanXchromosome>.4,"female",
#				ifelse(meanXchromosome<.38,"male","Unsure"))
#			datout=data.frame(datout,predictedGender=predictedGender,meanXchromosome=meanXchromosome)
			
#		} # end of if 
		
#	} # end of if

	datout
}



generate.aar <- function(agepred, phen)
{
	phen$index <- 1:nrow(phen)
	agepred <- subset(agepred, select=c(SampleID, DNAmAge))
	phen <- merge(phen, agepred, by.x="IID", by.y="SampleID", all.x=TRUE)

	message("The correlation between predicted age and actual age is ", cor(phen$Age_numeric, phen$DNAmAge, use="pair"))

	phen$AAR <- residuals(lm(Age_numeric ~ DNAmAge, phen, na.action=na.exclude))
	phen <- phen[order(phen$index), ]
	phen <- subset(phen, select=-c(index))
	return(phen)
}





##





# The function BMIQcalibration was created by Steve Horvath by heavily recycling code
# from A. Teschendorff's BMIQ function.
# BMIQ stands for beta mixture quantile normalization.
# Explanation: datM is a data frame with Illumina beta values (rows are samples, colums are CpGs.
# goldstandard is a numeric vector with beta values that is used as gold standard for calibrating the columns of datM.
# The length of goldstandard has to equal the number of columns of datM.
# Example code: First we impute missing values.
# library(WGCNA); dimnames1=dimnames(datMeth)
# datMeth= data.frame(t(impute.knn(as.matrix(t(datMeth)))$data))
# dimnames(datMeth)=dimnames1
# gold.mean=as.numeric(apply(datMeth,2,mean,na.rm=TRUE))
#datMethCalibrated=BMIQcalibration(datM=datMeth,goldstandard.beta=gold.mean)

BMIQcalibration=function(datM,goldstandard.beta,nL=3,doH=TRUE,nfit=20000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=FALSE,calibrateUnitInterval=TRUE){
	
	if (length(goldstandard.beta) !=dim(datM)[[2]] ) {stop("Error in function arguments length(goldstandard.beta) !=dim(datM)[[2]]. Consider transposing datM.")}
	if (plots ) {par(mfrow=c(2,2))}
	beta1.v = goldstandard.beta

	if (calibrateUnitInterval ) {datM=CalibrateUnitInterval(datM)}

	### estimate initial weight matrix from type1 distribution
	w0.m = matrix(0,nrow=length(beta1.v),ncol=nL);
	w0.m[which(beta1.v <= th1.v[1]),1] = 1;
	w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] = 1;
	w0.m[which(beta1.v > th1.v[2]),3] = 1;
	### fit type1
	print("Fitting EM beta mixture to goldstandard probes");
	set.seed(1)
	rand.idx = sample(1:length(beta1.v),min(c(nfit, length(beta1.v))  )   ,replace=FALSE)
	em1.o = blc(matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
	subsetclass1.v = apply(em1.o$w,1,which.max);
	subsetth1.v = c(mean(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]])),mean(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]])));
	class1.v = rep(2,length(beta1.v));
	class1.v[which(beta1.v < subsetth1.v[1])] = 1;
	class1.v[which(beta1.v > subsetth1.v[2])] = 3;
	nth1.v = subsetth1.v;
	print("Done");

	### generate plot from estimated mixture
	if(plots){
		print("Check");
		tmpL.v = as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
		tmpB.v = vector();
		for(l in 1:nL){
		  tmpB.v = c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
		}
		plot(density(beta1.v),main= paste("Type1fit-", sep=""));
		d.o = density(tmpB.v);
		points(d.o$x,d.o$y,col="green",type="l")
		legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
	}

	### Estimate Modes 
	if (  sum(class1.v==1)==1 ){ mod1U= beta1.v[class1.v==1]}
	if (  sum(class1.v==3)==1 ){ mod1M= beta1.v[class1.v==3]}
	if (  sum(class1.v==1) >1){ 
		d1U.o = density(beta1.v[class1.v==1])
		mod1U = d1U.o$x[which.max(d1U.o$y)]
	}
	if (  sum(class1.v==3)>1 ){ 
		d1M.o = density(beta1.v[class1.v==3])
		mod1M = d1M.o$x[which.max(d1M.o$y)]
	}

	### BETA 2
	for (ii in 1:dim(datM)[[1]] ){
		printFlush(paste("ii=",ii, "of", dim(datM)[[1]]))
		printFlush(date())
		sampleID=ii
		beta2.v = as.numeric(datM[ii,])

		d2U.o = density(beta2.v[which(beta2.v<0.4)]);
		d2M.o = density(beta2.v[which(beta2.v>0.6)]);
		mod2U = d2U.o$x[which.max(d2U.o$y)]
		mod2M = d2M.o$x[which.max(d2M.o$y)]

		### now deal with type2 fit
		th2.v = vector();
		th2.v[1] = nth1.v[1] + (mod2U-mod1U);
		th2.v[2] = nth1.v[2] + (mod2M-mod1M);

		### estimate initial weight matrix 
		w0.m = matrix(0,nrow=length(beta2.v),ncol=nL);
		w0.m[which(beta2.v <= th2.v[1]),1] = 1;
		w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] = 1;
		w0.m[which(beta2.v > th2.v[2]),3] = 1;

		print("Fitting EM beta mixture to input probes");
		# I fixed an error in the following line (replaced beta1 by beta2)
		set.seed(1)
		rand.idx = sample(1:length(beta2.v),min(c(nfit, length(beta2.v)),na.rm=TRUE)   ,replace=FALSE)
		em2.o = blc2(Y=matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol,verbose=TRUE);
		print("Done");

		### for type II probes assign to state (unmethylated, hemi or full methylation)
		subsetclass2.v = apply(em2.o$w,1,which.max);


		if (sum(subsetclass2.v==2)>0 ){
			subsetth2.v = c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),
			mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
		}
		if (sum(subsetclass2.v==2)==0 ){
			subsetth2.v = c(1/2*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 1/2*mean(beta2.v[rand.idx[subsetclass2.v==3]]), 1/3*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 2/3*mean(beta2.v[rand.idx[subsetclass2.v==3]]));
		}



		class2.v = rep(2,length(beta2.v));
		class2.v[which(beta2.v <= subsetth2.v[1])] = 1;
		class2.v[which(beta2.v >= subsetth2.v[2])] = 3;

		### generate plot
		if(plots){
			tmpL.v = as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
			tmpB.v = vector();
			for(lt in 1:nL){
				tmpB.v = c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
			}
			plot(density(beta2.v),  main= paste("Type2fit-",sampleID,sep="")  );
			d.o = density(tmpB.v);
			points(d.o$x,d.o$y,col="green",type="l")
			legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
		}

		classAV1.v = vector();classAV2.v = vector();
		for(l in 1:nL){
			classAV1.v[l] =  em1.o$mu[l,1];
			classAV2.v[l] =  em2.o$mu[l,1];
		}

		### start normalising input probes
		print("Start normalising input probes");
		nbeta2.v = beta2.v;
		### select U probes
		lt = 1;
		selU.idx = which(class2.v==lt);
		selUR.idx = selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
		selUL.idx = selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
		### find prob according to typeII distribution
		p.v = pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
		### find corresponding quantile in type I distribution
		q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
		nbeta2.v[selUR.idx] = q.v;
		p.v = pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
		### find corresponding quantile in type I distribution
		q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
		nbeta2.v[selUL.idx] = q.v;

		### select M probes
		lt = 3;
		selM.idx = which(class2.v==lt);
		selMR.idx = selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
		selML.idx = selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
		### find prob according to typeII distribution
		p.v = pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
		### find corresponding quantile in type I distribution
		q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
		nbeta2.v[selMR.idx] = q.v;


		if(doH){ ### if TRUE also correct type2 hemimethylated probes
			### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
			lt = 2;
			selH.idx = c(which(class2.v==lt),selML.idx);
			minH = min(beta2.v[selH.idx],na.rm=TRUE)
			maxH = max(beta2.v[selH.idx],na.rm=TRUE)
			deltaH = maxH - minH;
			#### need to do some patching
			deltaUH = -max(beta2.v[selU.idx],na.rm=TRUE) + min(beta2.v[selH.idx],na.rm=TRUE)
			deltaHM = -max(beta2.v[selH.idx],na.rm=TRUE) + min(beta2.v[selMR.idx],na.rm=TRUE)

			## new maximum of H probes should be
			nmaxH = min(nbeta2.v[selMR.idx],na.rm=TRUE) - deltaHM;
			## new minimum of H probes should be
			nminH = max(nbeta2.v[selU.idx],na.rm=TRUE) + deltaUH;
			ndeltaH = nmaxH - nminH;

			### perform conformal transformation (shift+dilation)
			## new_beta_H(i) = a + hf*(beta_H(i)-minH);
			hf = ndeltaH/deltaH ;
			### fix lower point first
			nbeta2.v[selH.idx] = nminH + hf*(beta2.v[selH.idx]-minH);

		}


		### generate final plot to check normalisation
		if(plots){
			print("Generating final plot");
			d1.o = density(beta1.v);
			d2.o = density(beta2.v);
			d2n.o = density(nbeta2.v);
			ymax = max(d2.o$y,d1.o$y,d2n.o$y);
			plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1), main=paste("CheckBMIQ-",sampleID,sep="") );
			points(d1.o$x,d1.o$y,col="red",type="l");
			points(d2n.o$x,d2n.o$y,col="blue",type="l");
			legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
		}

		datM[ii,]= nbeta2.v ;
	} # end of for (ii=1 loop
	datM
} # end of function BMIQcalibration





BMIQ = function(beta.v,design.v,nL=3,doH=TRUE,nfit=50000,th1.v=c(0.2,0.75),th2.v=NULL,niter=5,tol=0.001,plots=TRUE,sampleID=1,calibrateUnitInterval=TRUE){

	if (calibrateUnitInterval) {
		rangeBySample=range(beta.v,na.rm=TRUE)
		minBySample=rangeBySample[1]
		maxBySample=rangeBySample[2]
		if ( (minBySample<0 | maxBySample>1) & !is.na(minBySample) & !is.na(maxBySample) ) {
			y1=c(0.001,.999) 
			x1=c(minBySample,maxBySample)
			lm1=lm( y1 ~ x1 )
			intercept1=coef(lm1)[[1]]
			slope1=coef(lm1)[[2]]
			beta.v=intercept1+slope1*beta.v
		} # end of if
	} # end of if (calibrateUnitInterval


	type1.idx = which(design.v==1);
	type2.idx = which(design.v==2);

	beta1.v = beta.v[type1.idx];
	beta2.v = beta.v[type2.idx];


	### estimate initial weight matrix from type1 distribution
	w0.m = matrix(0,nrow=length(beta1.v),ncol=nL);
	w0.m[which(beta1.v <= th1.v[1]),1] = 1;
	w0.m[intersect(which(beta1.v > th1.v[1]),which(beta1.v <= th1.v[2])),2] = 1;
	w0.m[which(beta1.v > th1.v[2]),3] = 1;

	### fit type1
	print("Fitting EM beta mixture to goldstandard probes");
	set.seed(1)
	rand.idx = sample(1:length(beta1.v),min(c(nfit, length(beta1.v))  )   ,replace=FALSE)
	em1.o = blc2(Y=matrix(beta1.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
	subsetclass1.v = apply(em1.o$w,1,which.max);
	subsetth1.v = c(mean(max(beta1.v[rand.idx[subsetclass1.v==1]]),min(beta1.v[rand.idx[subsetclass1.v==2]])),mean(max(beta1.v[rand.idx[subsetclass1.v==2]]),min(beta1.v[rand.idx[subsetclass1.v==3]],na.rm=TRUE)));
	class1.v = rep(2,length(beta1.v));
	class1.v[which(beta1.v < subsetth1.v[1])] = 1;
	class1.v[which(beta1.v > subsetth1.v[2])] = 3;
	nth1.v = subsetth1.v;
	print("Done");

	### generate plot from estimated mixture
	if(plots){
		print("Check");
		tmpL.v = as.vector(rmultinom(1:nL,length(beta1.v),prob=em1.o$eta));
		tmpB.v = vector();
		for(l in 1:nL){
			tmpB.v = c(tmpB.v,rbeta(tmpL.v[l],em1.o$a[l,1],em1.o$b[l,1]));
		}

		pdf(paste("Type1fit-",sampleID,".pdf",sep=""),width=6,height=4);
		plot(density(beta1.v));
		d.o = density(tmpB.v);
		points(d.o$x,d.o$y,col="green",type="l")
		legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
		dev.off();
	}



	### Estimate Modes 
	if (  sum(class1.v==1)==1 ){ mod1U= beta1.v[class1.v==1]}
	if (  sum(class1.v==3)==1 ){ mod1M= beta1.v[class1.v==3]}
	if (  sum(class1.v==1) >1){ 
		d1U.o = density(beta1.v[class1.v==1])
		mod1U = d1U.o$x[which.max(d1U.o$y)]
	}
	if (  sum(class1.v==3)>1 ){ 
		d1M.o = density(beta1.v[class1.v==3])
		mod1M = d1M.o$x[which.max(d1M.o$y)]
	}


	d2U.o = density(beta2.v[which(beta2.v<0.4)]);
	d2M.o = density(beta2.v[which(beta2.v>0.6)]);
	mod2U = d2U.o$x[which.max(d2U.o$y)]
	mod2M = d2M.o$x[which.max(d2M.o$y)]


	### now deal with type2 fit
	th2.v = vector();
	th2.v[1] = nth1.v[1] + (mod2U-mod1U);
	th2.v[2] = nth1.v[2] + (mod2M-mod1M);

	### estimate initial weight matrix 
	w0.m = matrix(0,nrow=length(beta2.v),ncol=nL);
	w0.m[which(beta2.v <= th2.v[1]),1] = 1;
	w0.m[intersect(which(beta2.v > th2.v[1]),which(beta2.v <= th2.v[2])),2] = 1;
	w0.m[which(beta2.v > th2.v[2]),3] = 1;

	print("Fitting EM beta mixture to input probes");
	set.seed(1)
	rand.idx = sample(1:length(beta2.v),min(c(nfit, length(beta2.v)),na.rm=TRUE)   ,replace=FALSE)
	em2.o = blc2(Y=matrix(beta2.v[rand.idx],ncol=1),w=w0.m[rand.idx,],maxiter=niter,tol=tol);
	print("Done");

	### for type II probes assign to state (unmethylated, hemi or full methylation)
	subsetclass2.v = apply(em2.o$w,1,which.max);



	if (sum(subsetclass2.v==2)>0 ){
		subsetth2.v = c(mean(max(beta2.v[rand.idx[subsetclass2.v==1]]),min(beta2.v[rand.idx[subsetclass2.v==2]])),
		mean(max(beta2.v[rand.idx[subsetclass2.v==2]]),min(beta2.v[rand.idx[subsetclass2.v==3]])));
	}
	if (sum(subsetclass2.v==2)==0 ){
		subsetth2.v = c(1/2*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 1/2*mean(beta2.v[rand.idx[subsetclass2.v==3]]), 1/3*max(beta2.v[rand.idx[subsetclass2.v==1]])+ 2/3*mean(beta2.v[rand.idx[subsetclass2.v==3]]));
	}


	class2.v = rep(2,length(beta2.v));
	class2.v[which(beta2.v <= subsetth2.v[1])] = 1;
	class2.v[which(beta2.v >= subsetth2.v[2])] = 3;


	### generate plot
	if(plots){
		tmpL.v = as.vector(rmultinom(1:nL,length(beta2.v),prob=em2.o$eta));
		tmpB.v = vector();
		for(lt in 1:nL){
			tmpB.v = c(tmpB.v,rbeta(tmpL.v[lt],em2.o$a[lt,1],em2.o$b[lt,1]));
		}
		pdf(paste("Type2fit-",sampleID,".pdf",sep=""),width=6,height=4);
		plot(density(beta2.v));
		d.o = density(tmpB.v);
		points(d.o$x,d.o$y,col="green",type="l")
		legend(x=0.5,y=3,legend=c("obs","fit"),fill=c("black","green"),bty="n");
		dev.off();
	}

	classAV1.v = vector();classAV2.v = vector();
	for(l in 1:nL){
		classAV1.v[l] =  em1.o$mu[l,1];
		classAV2.v[l] =  em2.o$mu[l,1];
	}

	### start normalising input probes
	print("Start normalising input probes");
	nbeta2.v = beta2.v;
	### select U probes
	lt = 1;
	selU.idx = which(class2.v==lt);
	selUR.idx = selU.idx[which(beta2.v[selU.idx] > classAV2.v[lt])];
	selUL.idx = selU.idx[which(beta2.v[selU.idx] < classAV2.v[lt])];
	### find prob according to typeII distribution
	p.v = pbeta(beta2.v[selUR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
	### find corresponding quantile in type I distribution
	q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
	nbeta2.v[selUR.idx] = q.v;
	p.v = pbeta(beta2.v[selUL.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=TRUE);
	### find corresponding quantile in type I distribution
	q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=TRUE);
	nbeta2.v[selUL.idx] = q.v;

	### select M probes
	lt = 3;
	selM.idx = which(class2.v==lt);
	selMR.idx = selM.idx[which(beta2.v[selM.idx] > classAV2.v[lt])];
	selML.idx = selM.idx[which(beta2.v[selM.idx] < classAV2.v[lt])];
	### find prob according to typeII distribution
	p.v = pbeta(beta2.v[selMR.idx],em2.o$a[lt,1],em2.o$b[lt,1],lower.tail=FALSE);
	### find corresponding quantile in type I distribution
	q.v = qbeta(p.v,em1.o$a[lt,1],em1.o$b[lt,1],lower.tail=FALSE);
	nbeta2.v[selMR.idx] = q.v;


	if(doH){ ### if TRUE also correct type2 hemimethylated probes
		### select H probes and include ML probes (left ML tail is not well described by a beta-distribution).
		lt = 2;
		selH.idx = c(which(class2.v==lt),selML.idx);
		minH = min(beta2.v[selH.idx],na.rm=TRUE)
		maxH = max(beta2.v[selH.idx],na.rm=TRUE)
		deltaH = maxH - minH;
		#### need to do some patching
		deltaUH = -max(beta2.v[selU.idx],na.rm=TRUE) + min(beta2.v[selH.idx],na.rm=TRUE)
		deltaHM = -max(beta2.v[selH.idx],na.rm=TRUE) + min(beta2.v[selMR.idx],na.rm=TRUE)

		## new maximum of H probes should be
		nmaxH = min(nbeta2.v[selMR.idx],na.rm=TRUE) - deltaHM;
		## new minimum of H probes should be
		nminH = max(nbeta2.v[selU.idx],na.rm=TRUE) + deltaUH;
		ndeltaH = nmaxH - nminH;

		### perform conformal transformation (shift+dilation)
		## new_beta_H(i) = a + hf*(beta_H(i)-minH);
		hf = ndeltaH/deltaH ;
		### fix lower point first
		nbeta2.v[selH.idx] = nminH + hf*(beta2.v[selH.idx]-minH);

	}

	pnbeta.v = beta.v;
	pnbeta.v[type1.idx] = beta1.v;
	pnbeta.v[type2.idx] = nbeta2.v;

	### generate final plot to check normalisation
	if(plots){
		print("Generating final plot");
		d1.o = density(beta1.v);
		d2.o = density(beta2.v);
		d2n.o = density(nbeta2.v);
		ymax = max(d2.o$y,d1.o$y,d2n.o$y);
		pdf(paste("CheckBMIQ-",sampleID,".pdf",sep=""),width=6,height=4)
		plot(density(beta2.v),type="l",ylim=c(0,ymax),xlim=c(0,1));
		points(d1.o$x,d1.o$y,col="red",type="l");
		points(d2n.o$x,d2n.o$y,col="blue",type="l");
		legend(x=0.5,y=ymax,legend=c("type1","type2","type2-BMIQ"),bty="n",fill=c("red","black","blue"));
		dev.off();
	}

	print(paste("Finished for sample ",sampleID,sep=""));

	return(list(nbeta=pnbeta.v,class1=class1.v,class2=class2.v,av1=classAV1.v,av2=classAV2.v,hf=hf,th1=nth1.v,th2=th2.v));
}



CheckBMIQ = function(beta.v,design.v,pnbeta.v){### pnbeta is BMIQ normalised profile

	type1.idx = which(design.v==1);
	type2.idx = which(design.v==2);

	beta1.v = beta.v[type1.idx];
	beta2.v = beta.v[type2.idx];
	pnbeta2.v = pnbeta.v[type2.idx];
  
} # end of function CheckBMIQ


CalibrateUnitInterval=function(datM,onlyIfOutside=TRUE){

	rangeBySample=data.frame(lapply(data.frame(t(datM)),range,na.rm=TRUE))
	minBySample=as.numeric(rangeBySample[1,])
	maxBySample=as.numeric(rangeBySample[2,])
	if (onlyIfOutside) { 
		indexSamples=which((minBySample<0 | maxBySample>1) & !is.na(minBySample) & !is.na(maxBySample)) 
	}
	if (!onlyIfOutside) { indexSamples=1:length(minBySample)}
	if ( length(indexSamples)>=1 ){
		for ( i in indexSamples) {
			y1=c(0.001,0.999) 
			x1=c(minBySample[i],maxBySample[i])
			lm1=lm( y1 ~ x1 )
			intercept1=coef(lm1)[[1]]
			slope1=coef(lm1)[[2]]
			datM[i,]=intercept1+slope1*datM[i,]
		} # end of for loop
	}
	datM
} #end of function for calibrating to [0,1]


betaEst2=function (y, w, weights) 
{
	yobs = !is.na(y)
	if (sum(yobs) <= 1) 
		return(c(1, 1))
	y = y[yobs]
	w = w[yobs]
	weights = weights[yobs]
	N = sum(weights * w)
	p = sum(weights * w * y)/N
	v = sum(weights * w * y * y)/N - p * p
	logab = log(c(p, 1 - p)) + log(pmax(1e-06, p * (1 - p)/v - 
		1))
	if (sum(yobs) == 2) 
		return(exp(logab))
	opt = try(optim(logab, betaObjf, ydata = y, wdata = w, weights = weights, 
		method = "Nelder-Mead",control=list(maxit=50) ), silent = TRUE)
	if (inherits(opt, "try-error")) 
		return(c(1, 1))
	exp(opt$par)
} # end of function betaEst



blc2=function (Y, w, maxiter = 25, tol = 1e-06, weights = NULL, verbose = TRUE) 
{
	Ymn = min(Y[Y > 0], na.rm = TRUE)
	Ymx = max(Y[Y < 1], na.rm = TRUE)
	Y = pmax(Y, Ymn/2)
	Y = pmin(Y, 1 - (1 - Ymx)/2)
	Yobs = !is.na(Y)
	J = dim(Y)[2]
	K = dim(w)[2]
	n = dim(w)[1]
	if (n != dim(Y)[1]) 
		stop("Dimensions of w and Y do not agree")
	if (is.null(weights)) 
		weights = rep(1, n)
	mu = a = b = matrix(Inf, K, J)
	crit = Inf
	for (i in 1:maxiter) {
		warn0 = options()$warn
		options(warn = -1)
		eta = apply(weights * w, 2, sum)/sum(weights)
		mu0 = mu
		for (k in 1:K) {
			for (j in 1:J) {
				ab = betaEst2(Y[, j], w[, k], weights)
				a[k, j] = ab[1]
				b[k, j] = ab[2]
				mu[k, j] = ab[1]/sum(ab)
			}
		}
		ww = array(0, dim = c(n, J, K))
		for (k in 1:K) {
			for (j in 1:J) {
				ww[Yobs[, j], j, k] = dbeta(Y[Yobs[, j], j], 
				  a[k, j], b[k, j], log = TRUE)
			}
		}
		options(warn = warn0)
		w = apply(ww, c(1, 3), sum, na.rm = TRUE)
		wmax = apply(w, 1, max)
		for (k in 1:K) w[, k] = w[, k] - wmax
		w = t(eta * t(exp(w)))
		like = apply(w, 1, sum)
		w = (1/like) * w
		llike = weights * (log(like) + wmax)
		crit = max(abs(mu - mu0))
		if (verbose) 
			print(crit)
		if (crit < tol) 
			break
	}
	return(list(a = a, b = b, eta = eta, mu = mu, w = w, llike = sum(llike)))
}





main()
