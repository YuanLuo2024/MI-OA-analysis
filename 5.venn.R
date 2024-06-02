#install.packages("venn")


library(venn)
diffFile="diff.txt"
wgcnaFile="module_red.txt"
setwd("C:\\12.venn")
geneList=list()


rt=read.table(diffFile, header=T, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
uniqGene=unique(geneNames)
geneList[["Difference"]]=uniqGene

rt=read.table(wgcnaFile, header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])
uniqGene=unique(geneNames)
geneList[["WGCNA"]]=uniqGene


mycol=c("blue2","red2")
pdf(file="venn.pdf", width=5, height=5)
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F,ilabels=F)
dev.off()

intersectGenes=Reduce(intersect,geneList)
write.table(file="interGenes.txt", intersectGenes, sep="\t", quote=F, col.names=F, row.names=F)

