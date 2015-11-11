arguments <- commandArgs(T)
score <- arguments[1]
pheno <- arguments[2]
heightplotfile<-arguments[3]


r<-read.table(score,header=T,stringsAsFactors=F)
ph<-read.table(pheno,header=T)
m<-match(r$IID,ph$IID)
ph<-ph[m,]

pdf(heightplotfile,height=6,width=6)
plot(r$SCORE,ph$Height,cex=0.7,main=paste0("cor=",cor.test(ph$Height,r$SCORE)$estimate))
null <- dev.off()

message("Correlation between actual and predicted height: ", cor(ph$Height, r$SCORE, use="pair"))
