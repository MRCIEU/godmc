suppressMessages(library(meffil))

arguments <- commandArgs(T)
beta_file <- arguments[1]
prop_var <- as.numeric(arguments[2])
phen_file <- arguments[3]
pc_out <- arguments[4]


message("Loading methylation data")
load(beta_file)


message("Extracting most variable probes")
norm.beta <- norm.beta[meffil.most.variable.cpgs(norm.beta, n=20000), ]

message("Calculating PCs")
pc <- prcomp(t(norm.beta))

message("Identifying PCs that cumulatively explain ", prop_var, " of variance")
cumvar <- cumsum(pc$sdev^2) / sum(pc$sdev^2)
n_pcs <- which(cumvar > prop_var)[1]
pc <- pc$x[,1:n_pcs]

if(phen_file != "NULL")
{
	message("Removing PCs associated with EWAS phenotypes")
	phen <- read.table(phen_file, he=T)
	rownames(phen) <- phen$IID
	phen <- subset(phen, IID %in% rownames(pc), select=-c(IID))
	pc1 <- pc[rownames(pc) %in% rownames(phen), ]
	phen <- phen[match(rownames(pc1), rownames(phen)), ]
	stopifnot(all(rownames(phen) == rownames(pc1)))


	l <- lapply(1:ncol(phen), function(i)
	{
		pvals <- coefficients(summary(lm(phen[,i] ~ pc1)))[-1,4]
		which(pvals < 0.05/ncol(phen))
	})
	l <- sort(unique(unlist(l)))

	message("Identified ", length(l), " PC(s) associated with phenotypes")

	if(length(l) > 0)
	{
		pc <- pc[,! 1:ncol(pc) %in% l]
	}
}

pc1 <- t(pc)
save(pc, file=paste0(pc_out, ".RData"))
write.table(pc1, file=paste0(pc_out, ".txt"), row=T, col=T, qu=F, sep="\t")
