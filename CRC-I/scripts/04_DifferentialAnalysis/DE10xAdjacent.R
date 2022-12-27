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

label10x <- readRDS("data/cell/adjcent10x.rds")

subclstr10x <- SplitObject(label10x, split.by = 'Sub_Cluster')

# prepare to get DEGs numbers for each subcluster
#degnum <- tibble(id=1:length(subclstr10x), cellset = NA, degs = NA)

for (i in 26:length(subclstr10x)) {
  
  res <- tryCatch(deg(subclstr10x[i], id1='IT',id2='II'),
           error = function(e)
             return(NULL))
  if (is.null(res)){
    next
  }
  print(res[[1]])
  print(paste0(i, '/', length(subclstr10x)))
 
   write.csv(subset(res[[2]], p_val_adj < 0.05), file = paste0('results/DEG/adjacent10x/10xTumor-ITmarkers-',res[[1]], '.csv'))
  write.csv(subset(res[[2]]), file = paste0('results/DEG/adjacent10x/ALL10xTumor-ITmarkers-',res[[1]], '.csv'))
  png(filename = paste0('figures/volcano/adjacent10x/10xTumor-', res[[1]], '.png'), width = 600, height = 600)
  plot(res[[3]] +
         theme(text = element_text(size = 14))
  )
  dev.off()
  
  png(filename = paste0('figures/heatmap/10x/10xAdjacent-IT-', res[[1]], '.png'), width = 600, height = 600)
  plot(res[[4]]+
         scale_fill_viridis_c()+
         theme(text = element_text(size = 20)))
  dev.off()
  
  data <- res[[2]]
  genelist <- data[,2] # vector of fold change
  names(genelist) = rownames(data)# vector of gene id
  genelist = sort(genelist, decreasing = TRUE)
  
  gogsea <- gseGO(
    geneList = genelist,
    keyType = 'SYMBOL',
    eps=0,
    OrgDb = org.Hs.eg.db
  )
  rs <- tryCatch(dotplot(gogsea, showCategory=10),
           error=function(e)
             return(NULL))
  if(is.null(rs)){
    next
  }
  png(filename = paste0('figures/GSEAdot/10x/10xAdjacent-IT-', res[[1]], '.png'), width = 600, height = 600)
  plot(dotplot(gogsea, showCategory=10) + ggtitle(paste0(res[[1]], ' IT Gene Ontology')))
  dev.off()
}

for (i in 1:length(subclstr10x)){
  submark <- FindMarkers(subclstr10x[[i]], ident.1 = 'IT')
  print(names(subclstr10x[i]))
  print(dim(submark))
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
