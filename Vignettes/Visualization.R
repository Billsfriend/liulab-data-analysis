# Seurat official vignette
# https://satijalab.org/seurat/articles/visualization_vignette.html

library(Seurat)
library(tidyverse)
data("pbmc3k.final")
pbmc3k.final$groups <- sample(c("group1", "group2"), size = ncol(pbmc3k.final), replace = TRUE)
features <- c("LYZ", "CCL5", "IL32", "PTPRCAP", "FCGR3A", "PF4")
pbmc3k.final

RidgePlot(subclstr10x[["hB01_PlasmaB-IgG"]], features = 'FCGR2B')
VlnPlot(subclstr10x[["hB01_PlasmaB-IgG"]], features = 'FCGR2B')

select10x <- subclstr10x[c("hB01_PlasmaB-IgG",
                           "hB04_FollicularB-MS4A1",
                           'hM03_cDC2-CD1C',
                           'hM09_Macro-PLTP',
                           'hM10_Macro-IL1B',
                           'hM11_Monolike-FCN1',
                           'hM12_TAM-C1QC',
                           'hM13_TAM-SPP1')]
slct10x <- reduce(select10x, merge)
VlnPlot(labelSmart, group.by = 'Sub_Cluster',
        features = 'FCGR2B') + theme(legend.position = "none")
VlnPlot(label10x,
        features = 'FCGR2B', 
        group.by = 'Sub_Cluster',
        split.by = 'ITgeno',
        split.plot = TRUE)

VlnPlot(labelSmart,
        features = 'FCGR2B', 
        group.by = 'Sub_Cluster',
        split.by = 'ITgeno',
        split.plot = TRUE)
# Visualize the number of cell counts per sample
label10x@meta.data %>% 
  ggplot(aes(x=Sub_Cluster, fill = ITgeno)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

DotPlot(label10x,
        features = 'FCGR2B', 
        #group.by = 'Sub_Cluster',
        split.by = 'ITgeno')
# Single cell heatmap of feature expression
DoHeatmap(subset(label10x, downsample = 100),
          features = 'FCGR2B',
          group.by = 'Sub_Cluster',
          slot = 'data',
          size = 3) + theme(legend.position = "none")

# visualize taqman data
library(data.table)
library(tidyverse)
taqman <- fread('data/20220504-crc.txt') # edit x and y colname in advance
taqman <- taqman[Call != 'Negative'] # delete negative column
head(taqman)
ggplot(taqman, aes(x=FAM,y=VIC))+
  geom_point(aes(color=Call))+
  theme_classic()+
  scale_color_discrete(name='genotype', label=c('TT','II','IT'))
  theme(text = element_text(size=24))

ggplot(taqman, aes(x=Call))+
  geom_bar(aes(fill=Call))+
  geom_text(aes(label = ..count..), stat='count', vjust = 1.5, colour = "white")+
  theme_classic()+
  theme(text = element_text(size=16))
