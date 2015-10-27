removeOutliers <- function(x, thresh, remove=FALSE)
{
	m <- mean(x, na.rm=T)
	s <- sd(x, na.rm=T)
	index <- x > m + thresh*s | x < m - thresh*s

	if(remove)
	{
		x <- x[!index]
	} else {
		x[index] <- NA
	}
	return(x)
}


remRec <- function(x, thresh, iterations)
{
	d <- array(0, iterations+1)
	d[1] <- sum(!is.na(x))
	for(i in 1:iterations)
	{
		x <- removeOutliers(x, thresh)
		d[i+1] <- sum(!is.na(x))
	}
	return(list(x=x, its=d))
}



removeOutliersFromData <- function(X, thresh, iterations)
{
	n <- ncol(X)
	for(i in 1:n)
	{
		cat(i, "\n")
		X[,i] <- remRec(X[,i], thresh, iterations)$x
	}
	return(X)
}



arguments <- commandArgs(T)
pcafile <- arguments[1]
pcasd <- as.numeric(arguments[2])
npc <- as.numeric(arguments[3])
outliers <- arguments[4]
pcaplotfile<-as.character(arguments[5])

pca <- read.table(pcafile)
#pca<-read.table("./processed_data/genetic_data/pca.eigenvec")
#pcaplotfile=as.character("./processed_data/genetic_data/pcaplot.pdf")
#genetic_outliers="./processed_data/genetic_data/genetic_outliers.txt"

pca2 <- removeOutliersFromData(pca[,-c(1:2)], pcasd, iterations=3)
index <- apply(pca2, 1, function(x) any(is.na(x)))

genetic_outliers <- pca[index,1:2]

write.table(genetic_outliers, file=outliers, row=F, col=F, qu=F)

#if(length(genetic_outliers) > 0)
#{
#	write.table(subset(pca, !V2 %in% genetic_outliers[,2]), file=pcafile, row=F, col=F, qu=F)	
#}

library(ggplot2)
thresh1a <- mean(pca[,3],na.rm=T) + pcasd*(sd(pca[,3],na.rm=T))
thresh1b <- mean(pca[,3],na.rm=T) - pcasd*(sd(pca[,3],na.rm=T))
pc<-2:npc
scores<-data.frame()
hline.data1<-data.frame()
hline.data2<-data.frame()

out<-data.frame()
for (i in 1:length(pc)){
PC<-pc[i]+2

thresh2a <- mean(pca[,PC],na.rm=T) + pcasd*(sd(pca[,PC],na.rm=T))
thresh2b <- mean(pca[,PC],na.rm=T) - pcasd*(sd(pca[,PC],na.rm=T))

hline1 <- data.frame(PCy=paste("PC",pc[i],sep=""),z = thresh2a)
hline2 <- data.frame(PCy=paste("PC",pc[i],sep=""),z = thresh2b)

d <- data.frame(iid=pca$V2,threshold1=thresh1a,threshold2=thresh1b,threshold3=thresh2a,threshold4=thresh2b,PC=paste("PC1vsPC",pc[i],sep=""),PC1="PC1",PCy=paste("PC",pc[i],sep=""),PC.scores.x=pca[,3],PC.scores.y=pca[,PC])
#d2<-d[which(d[,1]%in%genetic_outliers[,2]),]

scores<-rbind(scores,d)
hline.data1<-rbind(hline.data1,hline1)
hline.data2<-rbind(hline.data2,hline2)

#d2<-d[d$PCy==2,]
#w1<-which(d2$PC.scores.y<thresh2b)
#w2<-which(d2$PC.scores.y>thresh2a)
#w<-unique(c(w1,w2))

#out<-rbind(out,d2)
out<-rbind(out,d)
}



vline.data1 <- data.frame(PCy=unique(scores$PCy),z = thresh1a)
vline.data2 <- data.frame(PCy=unique(scores$PCy),z = thresh1b)

    pcaplot <- ggplot(scores,aes(x=PC.scores.x, y=PC.scores.y)) +
    geom_point() +
    labs(y="PC.scores",x="PC.scores") +
    
    geom_vline(aes(xintercept = threshold1),colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept = threshold3),colour="blue", linetype="dashed") +
    geom_hline(aes(yintercept = threshold4),colour="blue", linetype="dashed") +
    geom_vline(aes(xintercept = threshold2),colour="blue", linetype="dashed") +
    facet_wrap(PC1~PCy ,ncol=4) +
    theme(strip.text.x = element_text(size=6),strip.text.y = element_text(size=6)) +
    theme(axis.text = element_text(size = 6)) 
    ggsave(plot=pcaplot,filename=pcaplotfile,height=6,width=6)



