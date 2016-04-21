arguments <- commandArgs(T)
related <- as.character(arguments[1])


message("Checking R version")
currentr <- paste0(R.Version()['major'], ".", R.Version()['minor'])
ch <- compareVersion(currentr, "3.1")
if(ch == -1)
{
	stop("You are running R version ", currentr, ". Please upgrade to at least 3.1.0.")
}



message("Checking that all required packages are present")

pkglist <- c(
	"lattice",
	"ggplot2",
	"data.table",
	"MatrixEQTL",
	"parallel",
	"matrixStats",
	"plyr",
	"meffil",
	"EasyQC",
	"impute"
)

index <- pkglist %in% rownames(installed.packages())
if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
}

pkglist_related<-c("GenABEL","SNPRelate","GENESIS")

if (related=="yes")
{
message("Checking that all required packages are present for related samples")
index <- pkglist_related %in% rownames(installed.packages())

if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
}
}

l<-list.files("./resources/genetics",pattern="1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz")
if(length(l)==0) 
{
    stop("Before continuing, you need to download 1kg_phase3_eur_aut_polymorphic.recoded.nodup.frq.gz from the sftp \n")    
}

message("All required packages are installed and required files are downloaded \n\n")
