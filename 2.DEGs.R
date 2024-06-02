#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")
#install.packages("ggplot2")


library(limma)
library(pheatmap)
library(ggplot2)

inputFile="geneMatrix.txt"
conFile="s1.txt"
treatFile="s2.txt"
logFCfilter=1
adj.P.Val.Filter=0.05

geoID="GSE"
conName="Control"
treatName="Treat"
setwd("C:\\05.diff\\GSE")

rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0,]

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
	rt[rt<0]=0
	rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum), rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut,file=paste0(geoID,".all.txt"),sep="\t",quote=F,col.names=F)

Type=c(rep(conName,conNum),rep(treatName,treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0(geoID,".diff.txt"),sep="\t",quote=F,col.names=F)

geneNum=50
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
    hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
    hmGene=diffGeneName
}
hmExp=data[hmGene,]
Type=c(rep(conName,conNum),rep(treatName,treatNum))
Type=factor(Type, levels=c(conName, treatName))
names(Type)=colnames(data)
Type=as.data.frame(Type)

pdf(file=paste0(geoID,".heatmap.pdf"), width=9, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         fontsize = 7,
         fontsize_row=5,
         fontsize_col=7)
dev.off()


allDiff$logFC[allDiff$logFC>20]=20
allDiff$logFC[allDiff$logFC< -20]=-20
Significant=ifelse((allDiff$adj.P.Val<adj.P.Val.Filter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")

p = ggplot(allDiff, aes(logFC, -log10(adj.P.Val)))+
    geom_point(aes(col=Significant))+
    scale_color_manual(values=c("green", "black", "red"))+
    labs(title = " ")+
    theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p=p+theme_bw()

pdf(paste0(geoID,".vol.pdf"),width=5.5,height=5)
print(p)
dev.off()

