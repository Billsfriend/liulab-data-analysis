library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)

data <- read.csv('results/DEG/tumor10x/10xTumor-ITmarkers-hT04_CD4-TCF7.csv')
data$entrez <- bitr(geneID = data[,1], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db, drop = FALSE)
data = subset(data, !is.na(data$entrez$ENTREZID))
genelist <- data[,3] # vector of fold change
names(genelist) = data$entrez$ENTREZID # vector of gene id
genelist = sort(genelist, decreasing = TRUE)

gogsea <- gseGO(
  geneList = genelist,
  keyType = 'ENTREZID',
  eps=0,
  pvalueCutoff = 0.05,
  OrgDb = org.Hs.eg.db,
)

dotplot(gogsea, showCategory=6) + ggtitle("GO for GSEA")

keggGsea <- gseKEGG(
  geneList = genelist,
  keyType = 'kegg',
  eps = 0,
  minGSSize = 3,
  scoreType = "pos"
)

dotplot(keggGsea, showCategory=6) + ggtitle("KEGG for GSEA")
