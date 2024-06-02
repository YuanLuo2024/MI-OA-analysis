##ggalluvial
library(ggalluvial)
setwd("C:\\Users\\ly\\Desktop\\sanky")
link <- read.table('input.txt', header = T, sep="\t", check.names=F)
x1 <- rep(link$TF[1:9],times=1092)  
x2 <- rep(link$TF[14:24],times=1092) 
x3 <- rep(link$TF[28:31],times=1092)
x <- c(x1,x2,x3)
y1 <- rep(link$`Hub Gene`[1],times=9828)  
y2 <- rep(link$`Hub Gene`[14],times=10920) 
y3 <- rep(link$`Hub Gene`[28],times=3276) 
y <- c(y1,y2,y3)
z1 <- rep(link$`Immune Cell`[1:13],times=756)  
z2 <- rep(link$`Immune Cell`[14:27],times=780)  
z3 <- rep(link$`Immune Cell`[28:39],times=273) 
z <- c(z1,z2,z3)

link <- data.frame(TFs=x,GENE=y,Cell=z) 

corLodes=to_lodes_form(link, axes = 1:ncol(link), id = "Cohort")

mycol <- rep(c("#F9837B","#F5886C","#F18C5A","#EC9042","#E79419","#E09719","#DA9C19","#D49F19","#CCA319","#C4A619",
               "#BBA919","#B1AC19","#A8B019","#9CB219","#8FB519","#81B819","#70BA19","#59BC19","#30BE19","#19C043",
               "#19C25A","#19C36B","#19C47A","#19C587","#19C694","#19C6A0","#19C7AB","#19C6B6","#19C6C1","#19C5CA",
               "#19C4D4","#19C3DC","#19C0E4","#19BEEC","#19BAF2","#19B7F9","#19B3FD","#19AEFF","#58A9FF","#7AA4FF",
               "#939EFF","#A798FF","#B892FF","#FF79A3","#FF76AF","#FF73BA","#FF71C4","#FF70CE","#FC70D8","#F971E0",
               "#F474E8","#EE78F0","#E77BF7","#DE81FC","#D386FF","#C68CFF","#A953FF","#FF4374"),58)

p1 <- ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) +  
  #用aes.flow控制线调颜色，forward说明颜色和前面一致，backward说明与后面一致
  geom_flow(width = 2/10,aes.flow = "forward",show.legend = TRUE) +  #画流动图
  geom_stratum(alpha = .9,width = 2/10, # 格子宽度
               linetype=1,size=0.4, # 格子的边框线
               color = "#696969",inherit.aes=TRUE) + #画冲击图
  scale_fill_manual(values = mycol) +
  #size = 2代表基因名字大小
  geom_text(stat = "stratum", size = 2.8,color="black") +
  xlab("") + ylab("") + theme_bw() + 
  theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  ggtitle("") + guides(fill = FALSE)   
ggsave(p1,filename = 'Sankey2.pdf',width = 8,height = 5)




link <- data.frame(
  TFs= c('STAT6','TP53','CREM','PGR','ATF2','E2F1'),
  GENE=c(c('DUSP1'),c('DUSP1'),c('DUSP1'),c('DUSP1'),c('DUSP1'),c('DUSP1')),
  Immune_cells= c(c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper'),
                  c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper'),
                  c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper'),
                  c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper'),
                  c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper'),
                  c('T_cells_CD4_memory_resting','Macrophages_M0','Mast_cells_resting','B_cells_memory','Plasma_cells','NK_cells_resting','Dendritic_cells_activated','Eosinophils','NK_cells_activated','Monocytes','Mast_cells_activated','Neutrophils','T_cells_follicular_helper')))
