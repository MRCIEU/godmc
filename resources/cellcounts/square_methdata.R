arguments <- commandArgs(T)

methylationfile <- arguments[1]
outfile <- arguments[2]

load(methylationfile)
norm.beta <- norm.beta^2

save(norm.beta, file=outfile)
