library(clusterProfiler)
library(tidyverse)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(pheatmap)

# read all degs from edgeR 
data <- read_csv("results/hisat-deg-entrez.csv")

goi <- read_csv('ref/GeneOfInterest.csv')

# select genes within our genes of interest list
dataOfInt <- filter(data, symbol %in% goi$symbol)
goi <- filter(goi, symbol %in% dataOfInt$symbol)

dataOfInt <- arrange(dataOfInt,symbol)
goi <- arrange(goi, symbol)

# add protein name for interest
dataOfInt$protein <- goi$protein

# tidy longer for ggplot
dataOfInt <- pivot_longer(dataOfInt, 9:10, names_to = 'type', values_to = 'CPM')

# extract genes of immunoglobulin heavy chain
igData <- data %>% 
  filter(symbol %in% c('Cd40','Ha')) %>%
  pivot_longer(9:10, names_to = 'type', values_to = 'CPM')

# list of mhc-2 genes
igData <- data %>% 
  filter(symbol %in% c('Cd40','H2-M3','H2-Oa','H2-M5','H2-Q4','H2-Aa','H2-DMb2','H2-DMa','H2-Ob','H2-T24','H2-Q5','H2-T22','Cd74')) %>%
  pivot_longer(9:10, names_to = 'type', values_to = 'CPM')

# plot goi logFC
dataOfInt %>%
  ggplot()+
  geom_col(aes(x=protein, y=logFC, fill=padj))+
  coord_flip()+
  scale_fill_viridis_c(name = 'P value')+
  theme_classic()+
  theme(text = element_text(size=16))+
  ylab('logFC (positive means upregulation in IgG-B cells)')

# plot cpm of B-T cell interaction candidates
dataOfInt %>%
  ggplot()+
  geom_col(aes(x=protein, y=log(CPM), fill = type), position = 'dodge')+
  coord_flip()+
  theme(text = element_text(size=16))

# plot cpm of immunoglobulin genes
igData %>%
  ggplot()+
  geom_col(aes(x=symbol, y=log(CPM), fill = type), position = 'dodge')+
  coord_flip()+
  xlab('Genes')+
  theme(text = element_text(size=16))

# plot ig-genes logFC
igData %>%
  ggplot()+
  geom_col(aes(x=symbol, y=logFC, fill=padj))+
  coord_flip()+
  scale_fill_viridis_c(name = 'P value')+
  theme_classic()+
  theme(text = element_text(size=16))+
  ylab('logFC (positive means upregulation in IgG-B cells)')+
  xlab('Genes')

# volcano plot
data <- data %>%
  mutate(genes=(padj<0.01 & abs(logFC)>2.5)*(logFC)/abs(logFC))
data$genes <- as.factor(data$genes)

data %>%
  ggplot()+
  geom_point(aes(x=logFC, y=-log10(padj), color=genes), size=1)+
  geom_vline(xintercept = c(-2.5,2.5), linetype='dashed')+
  geom_hline(yintercept = 2, linetype='dashed')+
  annotate(geom = 'text', x = -9, y = 59, label='Ighv11-2')+
  annotate(geom = 'text', x = 13, y = 53, label='Gm10052')+
  annotate(geom = 'text', x = 8, y = 48, label='Igkv9-124')+
  annotate(geom = 'text', x = -6, y = 42, label='Ly9c2')+
  scale_color_discrete(labels=c('IgM-UP','NS','IgG-UP'))

# uniqueData <- distinct(data, symbol, .keep_all = TRUE)

# all DEG as universe
universe <- data$entrez

# construct genelist with LFC numeric and symbol as names
genelist <- data$logFC
names(genelist) <- data$entrez

# extract sigDEG that padj < 0.05
siglist <- filter(data, padj<0.05)

Gup <- siglist %>%
  filter(logFC>0) %>%
  pull(entrez)

Gdown <-  siglist %>%
  filter(logFC<0) %>%
  pull(entrez)

#geneid <- bitr(names(genelist), OrgDb = org.Mm.eg.db, fromType = 'ENSEMBL', toType = 'ENTREZID', drop = TRUE)

genelist = sort(genelist, decreasing = TRUE)

# GO over-representation analysis of up-regulated in IgG-B
goEnrichUp <- enrichGO(gene          = Gup,
                     universe      = universe,
                     keyType = 'ENTREZID',
                     OrgDb         = org.Mm.eg.db,
                     ont           = "BP",
                     # minGSSize = 3,
                     # pvalueCutoff  = 0.01,
                     # qvalueCutoff  = 0.05, # stricter cutoff than default
                     readable      = TRUE)

dotplot(goEnrichUp, showCategory=6) +
  ggtitle("Upregulated pathways in IgG-B cells") +
  theme(text = element_text(size=16))

goplot(goEnrichUp, showCategory = 3)

cnetplot(goEnrichUp,
         showCategory = 4,
         foldChange = genelist,
         circular = TRUE,
         colorEdge = TRUE,
         node_label = 'category') +
  theme(text = element_text(size=14))

heatplot(goEnrichUp,
         foldChange = genelist,
         showCategory = 9)+
  theme(text = element_text(size=14))

goEnrichCls <- pairwise_termsim(goEnrich)
treeplot(goEnrichCls, showCategory = 20)

keggEnrich <- enrichKEGG(gene = )

# GO-ORA of up-regulated in IgM-B
goEnrichDown <- enrichGO(gene          = Gdown,
                         universe      = universe,
                         keyType = 'ENTREZID',
                         OrgDb         = org.Mm.eg.db,
                         ont           = "BP",
                         # minGSSize = 3,
                         # pvalueCutoff  = 0.01,
                         # qvalueCutoff  = 0.05, # stricter cutoff than default
                         readable      = TRUE)

dotplot(goEnrichDown, showCategory=8) +
  ggtitle("Upregulated pathways GO ORA in IgM-B cells") +
  theme(text = element_text(size=16))

cnetplot(goEnrichDown, showCategory = 4, foldChange = genelist, circular = TRUE, colorEdge = TRUE, node_label = 'category') + theme(text = element_text(size=14))

# heat plot show genes cluster of pathway
heatplot(goEnrichDown, foldChange = genelist, showCategory = 9)+ theme(text = element_text(size=14))

# similarity matrix
goEnrSimMat <- pairwise_termsim(goEnrichDown)

treeplot(goEnrSimMat,
         showCategory = 10,
         offset = 10)+
  theme(text = element_text(size=14))

# kegg ORA of IgM-B
MKeggEnrich <- enrichKEGG(gene = Gdown,
                         organism = 'mmu',
                         universe = universe)

dotplot(MKeggEnrich)

# GO GSEA
gogsea <- gseGO(
  geneList = genelist,
  #keyType = 'ENSEMBL',
  minGSSize = 3,
  eps=0,
  #pvalueCutoff = 0.2,
  OrgDb = org.Mm.eg.db,
  nPermSimple = 10000
)

# dot plot of GSEA
dotplot(gogsea, showCategory=8) +
  ggtitle("IgG-B cell GO:BP GSEA comparing with IgM-B cell") +
  theme(text = element_text(size=16))

# network plot show pathway nodes and genes
cnetplot(gogsea,
         showCategory = 4,
         foldChange = genelist,
         circular = TRUE,
         colorEdge = TRUE,
         node_label = 'category',
         readable = TRUE)

heatplot(gogsea,
         foldChange = genelist,
         showCategory = 9,
         readable = TRUE)+
  theme(text = element_text(size=14))

gogseaMat <- pairwise_termsim(gogsea)

treeplot(gogseaMat)

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

