library(Seurat)
library(tidyverse)
library(Platypus)
library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
source('Scripts/functions/GEX_volcano.R')

deg <- function(twoset, id1, id2){
  cname = names(twoset)
  
  submark <- FindMarkers(twoset[[1]], ident.1 = id1, ident.2 = id2, fc.name = 'avg_logFC')
  
  vol <- GEX_volcano2(DEGs.input = submark,
                      input.type = 'findmarkers',
                      condition.1 = id1,
                      condition.2 = id2,
                      maximum.overlaps = 6)
  
  # htmp <- GEX_DEgenes(GEX = twoset[[1]] ,min.pct = .25, FindMarkers.out = submark ,return.plot = "heatmap", group1 = 'IT', group2 = 'II', up.genes = 5,down.genes = 5,logFC = F)
  
  genelist <- submark[,2] # vector of fold change
  names(genelist) = rownames(submark)# vector of gene id
  genelist = sort(genelist, decreasing = TRUE)
  
  gogsea <- gseGO(
    geneList = genelist,
    keyType = 'SYMBOL',
    eps=0,
    pvalueCutoff = 0.05,
    OrgDb = org.Hs.eg.db
  )
  
  result = list(cname, submark, vol, gogsea)
}