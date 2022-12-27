# based on Seurat official DE test vignette

library(Seurat)
library(tidyverse)
library(Platypus)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
source('Scripts/functions/GEX_volcano.R')
source('Scripts/functions/DEG.R')

label10x <- readRDS("data/cell/label10x.rds")

label10x <- ScaleData(label10x)
label10x <- SetIdent(label10x, value = 'ITgeno')
#IT.markers <- FindMarkers(label10x, ident.1 = 'IT', fc.name = 'avg_logFC')
#head(IT.markers)
subclstr10x <- SplitObject(label10x, split.by = 'Sub_Cluster')
#subclstr10x <- subclstr10x[c('hB01','hB02','hB03','hB04','hB05','hM02','hM03','hM04','hM12','hM13')]
# 'hB09' cluster has no IT cells
write_rds(subclstr10x, 'data/cell/subclstr10x.rds')


# DEG identifing with only tumor infiltrates
subclstr10x <- read_rds('data/cell/tumor10x.rds')
subclstr10x <- SplitObject(subclstr10x, split.by = 'Sub_Cluster')

# prepare to get DEGs numbers for each subcluster
#degnum <- tibble(id=1:length(subclstr10x), cellset = NA, degs = NA)

for (i in 1:length(subclstr10x)) {
  
  res <- tryCatch(deg(subclstr10x[i], id1='IT',id2='II'),
           error = function(e)
             return(NULL))
  if (is.null(res)){
    next
  }
  print(res[[1]])
  print(paste0(i, '/', length(subclstr10x)))
  write.csv(subset(res[[2]], p_val_adj < 0.05), file = paste0('results/DEG/tumor10x/10xTumor-ITmarkers-',res[[1]], '.csv'))
  
  png(filename = paste0('figures/volcano/tumor10x/10xTumor-', res[[1]], '.png'), width = 600, height = 600)
  plot(res[[3]] +
         theme(text = element_text(size = 14))+
         annotate("segment",
                  x = 0.5, xend = 1,
                  y = 1, yend = 1,
                  arrow = arrow(ends = "last",
                                type = 'closed',
                                angle = 30,
                                length = unit(0.5,"cm")
                  )
         )+
         annotate('text',
                  x=0.9,y=2,label='IT',size=12) # add arrow and text indicating DEGs enriched in IT
  )
  dev.off()
  
  # png(filename = paste0('figures/heatmap/10x/10xTumor-IT-', res[[1]], '.png'), width = 600, height = 600)
  # plot(res[[4]]+scale_fill_viridis_c()+
  #        theme(text = element_text(size = 20)))
  # dev.off()
  
  data <- res[[2]]
  genelist <- data[,2] # vector of fold change
  names(genelist) = rownames(data)# vector of gene id
  genelist = sort(genelist, decreasing = TRUE)
  
  gogsea <- gseGO(
    geneList = genelist,
    keyType = 'SYMBOL',
    eps=0,
    pvalueCutoff = 0.05,
    OrgDb = org.Hs.eg.db
  )
  rs <- tryCatch(dotplot(gogsea, showCategory=10),
           error=function(e)
             return(NULL))
  if(is.null(rs)){
    next
  }
  png(filename = paste0('figures/GSEAdot/10x/10xTumor-IT-', res[[1]], '.png'), width = 600, height = 600)
  plot(dotplot(gogsea, showCategory=10) + ggtitle(paste0(res[[1]], ' IT Gene Ontology')))
  dev.off()
  
}

# bar plot: number of DEGs from each subcluster pairs
degnum <- subset(degnum, !is.na(degs) & degs != 0) # get rid of pairs with no significant DEGs

# define major clusters
degnum$bigset <- NA
degnum$bigset <- str_trunc(degnum$cellset, 2, 'right', ellipsis = '')

setDT(degnum)

ggplot(degnum)+
  geom_col(aes(x=cellset, y=degs, fill=bigset))+
  coord_flip()+
  ylab('Number of significant DEGs between II and IT')+
  scale_fill_discrete(
    name = 'Major clusters',
    limits = c("hT", "hM", "hI", 'hB'),
    labels = c("T cell", "Myeloid cell", "ILC", 'B cell')
  )

ggplot(degnum[bigset=='hT'])+
  geom_col(aes(x=cellset, y=degs, fill=bigset))+
  coord_flip()+
  ylab('Number of significant DEGs between II and IT')+
  scale_fill_discrete(
    name = 'Major clusters',
    labels = c('T cell')
  )+
  theme(text = element_text(size=16))
