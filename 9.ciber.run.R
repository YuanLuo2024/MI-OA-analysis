
#install.packages('e1071')
library('parallel')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("preprocessCore")

setwd("")
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "normalize.txt", perm=1000, QN=TRUE)

outTab=results[results[,"P-value"]<0.05,]
outTab=as.matrix(outTab[,1:(ncol(outTab)-3)])
outTab=rbind(id=colnames(outTab),outTab)
write.table(outTab, file="CIBERSORT-Results.txt", sep="\t", quote=F, col.names=F)