library("ggplot2")

data <- read.csv("expected.tsv", header=TRUE, sep="\t")
myPCA <- prcomp(data[, -1])
#plot(myPCA, type = "l")

#basic plot with eigenvectors
#biplot(myPCA, scale=0)

#extract PCA scores
#str(myPCA)


myPCA2 <- cbind(data, myPCA$x[,1:2])
pdf("testPca.pdf")
ggplot(myPCA2, aes(PC1, PC2, col = Sample, fill = Sample)) + geom_point(shape = 21, col = "black")+theme_classic()
dev.off()
