library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)

dir = c('C1/','C2/','OA1/','OA2/','OA3/','OA4/')
names(dir) = c('HC1', 'HC2','OA1','OA2','OA3','OA4')  
dir = c('C1/','C2/')
names(dir) = c('C1','C2/')
counts <- Read10X(data.dir =dir)  
counts <- Read10X('C1/')  
scRNA  <-  CreateSeuratObject(counts,min.cells = 3, min.features = 200)
scRNA$group <- scRNA$orig.ident
scRNA$group <- recode(scRNA$group,
                      HC1='HC',
                      HC2='HC',
                      OA1='OA',
                      OA2='OA',
                      OA3='OA',
                      OA4='OA')

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 

VlnPlot(scRNA,
              features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
              cols =brewer.pal(6,'Set3'),  
              pt.size = 0, 
              ncol = 2) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
scRNA1 <- subset(scRNA, subset = nFeature_RNA > 200& nFeature_RNA < 5000 & percent.mt < 10 & nCount_RNA < 50000)

library(harmony)
scRNA2 <- NormalizeData(scRNA1) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
system.time({scRNA_harmony <- RunHarmony(scRNA2, group.by.vars = "orig.ident")})

#scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.5)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30)
Idents(scRNA_harmony) <- 'seurat_clusters'
p1 <- DimPlot(scRNA_harmony, reduction = "umap",label = F) 
ggsave(p1,filename = 'OA-clusters.pdf',width = 8,height = 5)

library(SingleR)
hpca.se<-celldex::HumanPrimaryCellAtlasData() 
testdata <- GetAssayData(scRNA_harmony, slot="data")
clusters <- scRNA_harmony@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata,
                    ref = hpca.se, labels = hpca.se$label.main) 
celltype = data.frame(ClusterID=scRNA_harmony@meta.data$seurat_clusters, celltype=cellpred$labels, stringsAsFactors = F) 
scRNA_harmony@meta.data$celltype.singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
View(scRNA_harmony@meta.data)

a <- brewer.pal(10,'Paired')
DimPlot(scRNA_harmony, group.by="celltype.singleR", label=T, label.size=3, reduction="umap",split.by = 'group',cols = a)

#查看结果
Idents(scRNA_harmony) <- 'celltype.singleR'
DimPlot(scRNA_harmony)
DimPlot(scRNA_harmony,split.by = "group")
DimPlot(scRNA_harmony,group.by = "group")
p1=DimPlot(scRNA_harmony, group.by="seurat_clusters", split.by='group',label=T, label.size=4, reduction="umap")
p2=DimPlot(scRNA_harmony, group.by="celltype.singleR", split.by='group',label=T, label.size=4, reduction="umap")
p3=DimPlot(scRNA_harmony,group.by = "group",label = T)
plota=p1+p2+p3
plota
ggsave("singleR.pdf", plot=plota, width=12, height=6) 
View(scRNA_harmony@meta.data)

Idents(scRNA_harmony) <- 'seurat_clusters'
markers <- FindAllMarkers(object = scRNA_harmony, 
                          only.pos = TRUE,
                          logfc.threshold = 0.25) 
save(markers,file = 'markers.Rdata')
DefaultAssay(scRNA_harmony) <- "RNA" 

all.markers = markers %>% select(gene, everything()) %>% subset(p_val_adj<0.05)
top10.markers = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

scRNA_harmony$celltype.cellmarker <- scRNA_harmony$seurat_clusters
scRNA_harmony$celltype.cellmarker <- recode(scRNA_harmony$celltype.cellmarker,
                                            '0' = 'Chondrocytes',
                                            '1' = 'Chondrocytes',
                                            '2' = 'Chondrocytes',
                                            '3' = 'Chondrocytes',
                                            '4' = 'Chondrocytes',
                                            '5' = 'Chondrocytes',
                                            '6' = 'Chondrocytes',
                                            '7' = 'Chondrocytes',
                                            '8' = 'Chondrocytes',
                                            '9' = 'Chondrocytes',
                                            '10' = 'Chondrocytes',
                                            '11' = 'Chondrocytes',
                                            '12' = 'Endothelial_cells',
                                            '13' = 'Endothelial_cells',
                                            '14' = 'B_cells')

Idents(scRNA_harmony) <- 'seurat_clusters' 
gene <- c('SOX9','COL2A1','VWF','SPARCL1','PTPRC','CXCR4')
p1 <- DotPlot(scRNA3, features = gene ,
        assay='RNA' ) + 
  coord_flip() + 
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(angle = 45, hjust = 0.5,vjust=0.5))+ 
  labs(x=NULL,y=NULL) + 
  guides(size = guide_legend("Percent Expression") )+ #legend
  scale_color_gradientn(colours = c('#330066','#336699','#66CC66','#FFCC33')) 
ggsave('markers.pdf',p1,width = 8,height = 6)

a <- c('DUSP1', 'FOS', 'THBS1')
dents(scRNA3) <- 'group'
scRNA_OA <- subset(scRNA_harmony,idents = 'OA')
scRNA_HC <- subset(scRNA_harmony,idents = 'HC')

for(i in 1:length(a)){
  p1 <- plot_density(scRNA_OA, a[i])
  ggsave(filename = paste(a[i],'-density(OA).pdf',sep = ''),plot = p1, width = 6,height = 4 )
} 

for(i in 1:length(a)){
  p1 <- plot_density(scRNA_MI, a[i])
  ggsave(filename = paste(a[i],'-density(MI).pdf',sep = ''),plot = p1, width = 6,height = 4 )
} 

cols = c('#67a7bd','#edbac8')
my_comparisons <- list(c('HC','OA'))
for(i in 1:length(a)){
  scRNA_harmony <- AddModuleScore(scRNA_harmony,features = a[i],name = paste(a[i],'_score',sep = ''))
  p1 <- ggviolin(scRNA_harmony@meta.data, x = "group", y = paste(a[i],'_score1',sep = ''),
                 color = "group",add = 'mean_sd',fill = 'group',
                 add.params = list(color = "black")) + 
    stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
    scale_color_manual(values = cols) + 
    scale_fill_manual(values =  cols) +
    theme(axis.text.x.bottom = element_text(angle = 0,vjust = 0.5,hjust = 1)) + 
    NoLegend() + labs(x = '')
  ggsave(filename = paste(a[i],'-violin.pdf',sep = ''),plot = p1, width = 6,height = 4 )
} 

