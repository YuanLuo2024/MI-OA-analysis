library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(RColorBrewer)
library(ggpubr)
library(Nebulosa)

load('scRNA.Rdata')

Idents(scRNA_harmony) <- 'cell_type_original' 
gene <- c('PLIN1','PNPLA2','GPAM',
          'TNNT2','MYO18B','MYBPC3',
          'TOP2A',
          'VWF','LDB2','EGFL7',
          'PDGFRA','FBN1','COL5A1',
          'IL7R','BCL11B',
          'KIT','SLC24A3','CPA3',
          'RBM47','CD14','CD163',
          'NRXN1','CADM2','CDH19',
          'KCNJ8','NOTCH3',
          'MYH11','MYL9','ACTA2')
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
Idents(scRNA3) <- 'major_labl'
scRNA_MI <- subset(scRNA3,idents = 'IZ')
scRNA_HC <- subset(scRNA3,idents = 'CTRL')

for(i in 1:length(a)){
  p1 <- plot_density(scRNA_MI, a[i])
  ggsave(filename = paste(a[i],'-density(MI).pdf',sep = ''),plot = p1, width = 6,height = 4 )
} 

for(i in 1:length(a)){
  p1 <- plot_density(scRNA_HC, a[i])
  ggsave(filename = paste(a[i],'-density(HC).pdf',sep = ''),plot = p1, width = 6,height = 4 )
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

