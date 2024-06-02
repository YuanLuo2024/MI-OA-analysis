#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")
#install.packages("circlize")
#install.packages("RColorBrewer")
#install.packages("ggpubr")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("ComplexHeatmap")


library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(ComplexHeatmap)

pvalueFilter=0.05
qvalueFilter=1

colorSel="qvalue"
if(qvalueFilter>0.05){
	colorSel="pvalue"
}

setwd("C:\\07.GO")
rt=read.table("interGenes.List.txt", header=F, sep="\t", check.names=F)

genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
gene=entrezIDs[entrezIDs!="NA"]
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)

kk=enrichGO(gene=gene, OrgDb=org.Hs.eg.db, pvalueCutoff=1, qvalueCutoff=1, ont="all", readable=T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO, file="GO.txt", sep="\t", quote=F, row.names = F)

pdf(file="barplot.pdf", width=9, height=7)
bar=barplot(kk, drop=TRUE, showCategory=10, label_format=100, split="ONTOLOGY", color=colorSel) + facet_grid(ONTOLOGY~., scale='free')
print(bar)
dev.off()
		