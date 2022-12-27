library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(pheatmap)

# filter only DEGs with adj.p<0.05 as enriched genes 
data <- read.csv("DiffAnalysis/DEG_SalTx_0.05.csv", row.names = 1)

# all DEG as universe
universe <- read.csv("DiffAnalysis/DEG_SalTx_results.all.csv", row.names = 1)

genelist <- data[,2]
unilist <- universe[,2]
names(genelist) <- rownames(data)
#geneid <- bitr(names(genelist), OrgDb = org.Mm.eg.db, fromType = 'ENSEMBL', toType = 'ENTREZID', drop = TRUE)

genelist = sort(genelist, decreasing = TRUE)
unilist = sort(unilist, decreasing = TRUE)

over75 <- names(genelist)[abs(genelist) > 0.75]
over75p <- names(genelist)[genelist > 0.75]
over75m <- names(genelist)[genelist < -0.75]
over200 <- names(genelist)[abs(genelist) > 2]

# GO over-representation analysis
goEnrich <- enrichGO(gene          = over75,
                universe      = names(unilist),
                keyType = 'ENSEMBL',
                OrgDb         = org.Mm.eg.db,
                ont           = "BP",
                # pvalueCutoff  = 0.01,
                # qvalueCutoff  = 0.05, # stricter cutoff than default
                readable      = TRUE)

dotplot(goEnrich, showCategory=9) + ggtitle("Differential pathways in Knock-out") + theme(text = element_text(size=16))

goplot(goEnrich, showCategory = 3)

cnetplot(goEnrich, showCategory = 4, foldChange = genelist, circular = TRUE, colorEdge = TRUE, node_label = 'category') + theme(text = element_text(size=14))

heatplot(goEnrich, foldChange = genelist, showCategory = 9)+ theme(text = element_text(size=14))

goEnrichCls <- pairwise_termsim(goEnrich)
treeplot(goEnrichCls, showCategory = 20)

keggEnrich <- enrichKEGG(gene = )

gogsea <- gseGO(
  geneList = genelist,
  keyType = 'ENSEMBL',
  minGSSize = 3,
  #eps=0,
  pvalueCutoff = 0.2,
  OrgDb = org.Mm.eg.db,
)

dotplot(gogsea, showCategory=6) + ggtitle("GO for GSEA")

dataup <- subset(data, data$log2FoldChange >0)
datadown <- subset(data, data$log2FoldChange <0)
updown <- merge(head(dataup), head(datadown), all=TRUE)
write.csv(updown, 'DiffAnalysis/updown.csv')
# edit in excel...
updownp <- read.csv('DiffAnalysis/heatmap.txt', sep = '\t')
updownp$symbol <- reorder(updownp$symbol, updownp$factor, mean)
# heatmap
ggplot(updownp, aes(x=sample, y=symbol, fill=relativeExp))+
  geom_tile()+
  geom_raster()+
  scale_fill_viridis_c(option = 'B')+
  theme(text = element_text(size=16))


# volcano
alldata <- read.csv('DiffAnalysis/DEG_SalTx_results.all.csv')
alldata$minusLogQ <- -log(alldata$padj)
alldata$highlight <- NA
alldata$highlight[which(alldata$padj<0.001 & alldata$log2FoldChange > 1)] <- 'up'
alldata$highlight[which(alldata$padj<0.001 & alldata$log2FoldChange < -1)] <- 'down'
ggplot(alldata, aes(x=log2FoldChange, y=minusLogQ, color=highlight))+
  geom_point()+
  ylab('-logQvalue')+
  scale_color_discrete(labels = c('downregulated in CKO','upregulated in CKO', 'NS'))+
  annotate("text", x = 5.5, y = 250, label = "Hspg2")+ 
  annotate("text", x = 6.3, y = 150, label = "Vegfa")+
  annotate("text", x = -3.5, y = 110, label = "Cd79a")+ # manually label high point
  theme(text = element_text(size = 16))
