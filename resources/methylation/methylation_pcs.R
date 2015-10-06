library(meffil)

arguments <- commandArgs(T)
beta_file <- arguments[1]
n_pcs <- as.numeric(arguments[2])
pc_out <- arguments[3]


load(beta_file)

norm.beta <- norm.beta[meffil.most.variable.cpgs(norm.beta, n=20000), ]
pc <- prcomp(t(norm.beta))$x[,1:n_pcs]
pc1 <- t(pc)
save(pc, file=paste0(pc_out, ".RData"))
write.table(pc1, file=paste0(pc_out, ".txt"), row=T, col=T, qu=F, sep="\t")