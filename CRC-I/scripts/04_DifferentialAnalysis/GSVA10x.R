# ref: https://www.jianshu.com/p/adda4536b2cb

library(GSVA)
library(GSVAdata)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(msigdbr)
library(limma) #limma做差异分析

data(c2BroadSets)

# subclstr10x <- readRDS("data/cell/subclstr10x.rds")
# try only analyze tumor infiltrates here
tumor10x <- read_rds('data/cell/tumor10x.rds')
tumor10x <- SetIdent(tumor10x, value = "ITgeno")
subclstr10x <- SplitObject(tumor10x, split.by = 'Sub_Cluster')

# extract only cell subset of interest (high expression of FcgR2b)
fcgrset <- read_csv('results/fcgr2bexpr.txt', col_names = 0)

fcgr10x <- subclstr10x[fcgrset[[1]]]
write_rds(fcgr10x,'data/cell/fcgr10x.rds')

m_df = msigdbr(species = "Homo sapiens",
               category = "C2",
               subcategory = 'CP:REACTOME') # C3, C5-GO are worth trying
#m_df = subset(m_df, gs_subcat != 'CGP')
# C2 (-CGP) is still too big, just CP:reactome?
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

# mygsva(fcgr10x,msigdbr_list)
# 
# mygsva <- function(seulist, msigdbr_list){
#   for (i in 1:length(seulist)){
#     scl <- seulist[[i]]
#     cname <- names(seulist[i])
#     print(paste0(i, '/', length(seulist),': ',cname))
#     # compare 2 groups
#     DD1 <-  subset(scl, idents = c("II", "IT")) # control group first
#     expr <- as.data.frame(DD1@assays$RNA@data)
#     meta <- DD1@meta.data[, c("ITgeno")]
#     expr = as.matrix(expr)
#     kegg <- gsva(expr, msigdbr_list,
#                  parallel.sz = 10, min.sz = 10, max.sz = 100)
#     group <-
#       factor(meta, levels = c("II", "IT"), ordered = F)## put control group first
#     design <- model.matrix( ~ group)# 构建比较矩阵
#     colnames(design) <- levels(group)
#     fit <- lmFit(kegg, design)#线性模型拟合
#     fit1 <- eBayes(fit)#贝叶斯检验
#     allDiff = topTable(
#       fit1,
#       adjust = 'fdr',
#       coef = 2,
#       number = Inf,
#       resort.by = 't',
#       p.value = 0.05
#     )
#     write.csv(allDiff, paste0("results/GSVA/tumor10x",'-', cname, '.csv'))
#   }
# }

# scripts below only generate strict GSVA csv
# subcluster 27-hM09 have no IT, skip
# in tumor infiltrates, cluster 6,14,20,37,38 have low cell 
for (i in 12:length(fcgr10x)){
  scl <- fcgr10x[[i]]
  cname <- names(fcgr10x[i])
  print(paste0(i, '/', length(fcgr10x),': ',cname))
  # compare 2 groups
  DD1 <-  subset(scl, idents = c("II", "IT")) # control group first
  expr <- as.data.frame(DD1@assays$RNA@data)
  meta <- DD1@meta.data[, c("ITgeno")]
  expr = as.matrix(expr)
  kegg <- gsva(expr, msigdbr_list,
               parallel.sz = 10, min.sz = 10, max.sz = 100)
  group <-
    factor(meta, levels = c("II", "IT"), ordered = F)## put control group first
  design <- model.matrix( ~ group)# 构建比较矩阵
  colnames(design) <- levels(group)
  fit <- lmFit(kegg, design)#线性模型拟合
  fit1 <- eBayes(fit)#贝叶斯检验
  allDiff = topTable(
    fit1,
    adjust = 'fdr',
    coef = 2,
    number = Inf,
    resort.by = 't',
    p.value = 0.05 # or 0.01 for a stringent one
  )
  write.csv(allDiff, paste0("results/GSVA/cell/Stct10x-", cname, '.csv'))
}

# plot digested pathways for GSVA
library(data.table)
library(tidyverse)

richpaths <- fread('results/GSVA/digest/IIITall.txt')

richpath <- richpaths[subset=='CXCR6+ CD4+ RM']

# reorder by t value
richpath$pathway <- reorder(richpath$pathway, richpath$t, mean)

# draw barplot of GSVA t value
ggplot(richpath, aes(pathway, t))+
  geom_col(aes(fill=factor)) +
  theme_classic()+
  theme(
    text = element_text(size = 16), 
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  geom_text(aes(label=pathway, hjust=0.5+0.5*t/abs(t)), size=5,nudge_y = -1*(richpath$t)-0.1*richpath$t/abs(richpath$t))+
  coord_flip() +
  geom_hline(yintercept=0)+
  ylab('t value of GSVA analysis')+
  ggtitle(richpath$subset[1])

# proportion of T cell subset in all T cells
natperc <- read_delim('results/GSVA/digest/naturepercent.txt')

natperc$subset <- reorder(natperc$subset, natperc$mean, mean)

nperc <- ggplot(natperc,aes(subset,mean))
nperc+
  geom_col(aes(fill=factor))+
  coord_flip()

# generate heatmap by geom_tile
ggplot(natperc, aes(x=factor,y=subset,fill=mean))+
  geom_tile()+
  geom_raster()+
  scale_fill_viridis_c()+
  labs(fill='Proportion in all T cells')+
  theme(text=element_text(size=16))
