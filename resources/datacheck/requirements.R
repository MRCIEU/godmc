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
	"GenABEL",
	"matrixStats",
	"plyr",
	"SNPRelate",
	"GENESIS",
	"meffil"
)

index <- pkglist %in% rownames(installed.packages())
if(any(!index))
{
	stop("Before continuing, the following packages need to be installed:\n", paste(pkglist[!index], collapse="\n"))
}

message("All required packages installed\n\n")
