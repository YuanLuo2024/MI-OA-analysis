#install.packages("circlize")
#install.packages("RColorBrewer")


library(circlize)
library(RColorBrewer)

input="KEGG.txt"
outpdf="KEGG.pdf"
setwd("D:\\KEGG")


data=read.table(input, header=T, sep="\t", check.names=F)
data=data[order(data$pvalue),]
data=head(data, n=8)
pvalue=round(-log(data$pvalue,10),2)
genelist=unlist(strsplit(data$geneID,"/"))
genetable=table(genelist)

pathway_col = colorRampPalette(brewer.pal(9, "Paired"))(nrow(data))
genes_col = rainbow(length(names(genetable)),s=0.7,v=0.7)
pvalue_col_fun <- colorRamp2(
	breaks = c(min(pvalue), mean(pvalue), max(pvalue)), 
	colors = c("#FFFFCC", "#FF9966", "#FF0000")
)

genedict = list()
for(i in names(genetable)) genedict[[i]] = 0
pathwaydict = list()
for(i in data$Description) pathwaydict[[i]] = 0
# link  gene--pathway
gene_pathway_list = list()
n = 0
for(i in 1:nrow(data)){
  genei = strsplit(data$geneID[i],"/")[[1]] # multi
  pathwayi = data$Description[i]
  for(j in 1:length(genei)){
    n = n+1
    genej = genei[j] # single
    gene_pathway_list[[n]] = data.frame(gene=genej,pathway=data$Description[i],start=genedict[[genej]],
                                        end=genedict[[genej]]+1,start2 = pathwaydict[[pathwayi]],end2=pathwaydict[[pathwayi]]+1,
                                        pvalue=pvalue[i])
    genedict[[genej]] = genedict[[genej]]+1
    pathwaydict[[pathwayi]] = pathwaydict[[pathwayi]]+1
  }
}
#gene \t pathway \t start \t end \t pvalue
gene_pathway_data = as.data.frame(do.call('rbind',gene_pathway_list))
gene_pathway_data$linkcol = pathway_col[as.numeric(as.factor(gene_pathway_data$pathway))]
#right pathway
data3 = data.frame(id=data$Description,start = 0, end = data$Count)
# left top
data1 = data.frame(id=names(genetable),start=0,end = as.numeric(genetable))
# main chrom
df = as.data.frame(rbind(data3,data1))

get_sig = function(p){
  ifelse(p> -log(0.001,10),"***",ifelse(p> -log(0.01,10),'**','*'))
}

bed3 = data.frame(data3,yt=0,yb=1,col=pathway_col[as.numeric(as.factor(data3$id))],p=0,text='')
bed1 = data.frame(data1,yt=0.5,yb=1,col=genes_col[as.numeric(as.factor(data1$id))],p=0,text='')
bed2 = data.frame(id=gene_pathway_data$gene,start=gene_pathway_data$start,
                  end=gene_pathway_data$end,yt=0,yb=0.5,
                  col=pvalue_col_fun(gene_pathway_data$pvalue),p=gene_pathway_data$pvalue,
                  text=get_sig(gene_pathway_data$pvalue))

bed = as.data.frame(rbind(bed1,bed2,bed3))

pdf(file=outpdf, width=12, height=7)
layout(mat=matrix(c(1,1,1,0,2,3),nc=2),width=c(6.5,3.5),height=c(2,2,7))
circos.par(track.margin=c(0.01,0.01), start.degree=90)
circos.genomicInitialize(df,plotType="none")

circos.trackPlotRegion(ylim = c(0, 1), panel.fun = function(x, y) {
  sector.index = get.cell.meta.data("sector.index")
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  if(!any(data3$id%in%sector.index)){
    circos.text(mean(xlim), mean(ylim), sector.index, cex = 1, facing = "bending.inside", niceFacing = TRUE)
  }
}, track.height = 0.08, bg.border = NA,bg.col = NA)

circos.genomicTrack(bed, ylim = c(0, 1),track.height = 0.15,bg.border=NA,
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicRect(region, value, ytop = value[,1], ybottom = value[,2], col = value[,3],
                                         border = "black", ...)
                      for(j in 1:nrow(value)){
                        if(value[j,4]!=0){
                          circos.genomicText(region[j,], value[j,], y = 0.25, labels = value[j,5], adj=0.5,cex=1,...)
                        }
                      }
                    })


for(i in 1:nrow(gene_pathway_data)){
  genei = gene_pathway_data$gene[i]
  pathwayi = gene_pathway_data$pathway[i]
  circos.link(genei, c(gene_pathway_data$start[i], gene_pathway_data$end[i]), 
              pathwayi, c(gene_pathway_data$start2[i], gene_pathway_data$end2[i]), 
              col = gene_pathway_data$linkcol[i], border = "black")
}
circos.clear()


par(mar=c(3,0,3,10))
barplot(rep(1,100),col=pvalue_col_fun(seq(min(pvalue),max(pvalue),length=100)),
        space=0,border=NA,xaxt="n",yaxt="n",main="-log10(pvalue)")
axis(1,c(1,50,100),c(min(pvalue),mean(pvalue),max(pvalue)),tick=F)


par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,xlab="",ylab="")
legend("left",legend=bed3$id,col=bed3$col,pch=15,pt.cex=3,cex=1.2,bty="n",title="KEGG pathway")
dev.off()

