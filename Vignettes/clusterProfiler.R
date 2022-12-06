# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
# BiocManager::install('clusterProfiler')
# BiocManager::install(version = '3.15') for R 4.2
library(clusterProfiler)
library(tidyverse)
library(org.Hs.eg.db)

# GO classification
# use genelist in DOSE package as example
# contain 2 columns: gene ids and fold change number
data(geneList, package="DOSE")
folc <- as.data.frame(geneList)
# select genes with significant change (>2)
gene <- names(geneList)[abs(geneList) > 2]

# Entrez gene ID
head(gene)

ggo <- groupGO(gene     = gene,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE) # readable = T convert gene ID to symbol

head(ggo)

# GO over-representation analysis
ego <- enrichGO(gene          = gene,
                universe      = names(geneList),
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)
# Any gene ID type that is supported in OrgDb can be directly used in GO analyses. Users need to specify the keyType parameter to specify the input gene ID type.
head(ego)

# GO GSEA
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              eps = 0,
              verbose      = TRUE)

# Visualize enriched GO terms as a directed acyclic graph
goplot(ego)

# GO semantic similarity can be calculated by GOSemSim. 
# We can use it to cluster genes/proteins into different clusters based on their functional similarity 
# and can also use it to measure the similarities among GO terms to reduce the redundancy of GO enrichment results

# KEGG enrichment analysis
# select supported organism
ecoli <- search_kegg_organism('Homo sapiens', by='scientific_name')
# output contain KEGG_code, scientific name, and commmon name
head(ecoli)
kk <- enrichKEGG(gene         = gene,
                 organism     = 'hsa', # human as KEGG code
                 pvalueCutoff = 0.05)
head(kk)

# KEGG pathway GSEA
kk2 <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
               minGSSize    = 120,
               pvalueCutoff = 0.05,
               verbose      = TRUE)
head(kk2)

# KEGG module sometimes offer a more straightforward interpretation
mkk <- enrichMKEGG(gene = gene,
                   organism = 'hsa',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)                   

mkk2 <- gseMKEGG(geneList = geneList,
                 organism = 'hsa')
dim(mkk2)

# visualize enriched KEGG pathway
# BiocManager::install('pathview')
library(pathview)
hsa04110 <- pathview(gene.data  = geneList,
                     pathway.id = "hsa04110",
                     species    = "hsa",
                     limit      = list(gene=max(abs(geneList)), cpd=1),
                     kegg.dir = 'figures/kegg')
# return pathway plot png file in working directory

# wikipathways analysis
enrichWP(gene, organism = "Homo sapiens")
gseWP(geneList, organism = "Homo sapiens")
# no visualization?

# reactome analysis
# BiocManager::install('ReactomePA')

# bitr() function translate multiple names
eg = bitr(gene, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")

#
# visualization
library(DOSE)
library(enrichplot)

edo <- enrichDGN(gene) 
#Enrichment analysis based on the DisGeNET

barplot(edo, showCategory=10) 

# change y axis value to qscore
mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

# Dot plot is similar to bar plot with the capability to encode another score as dot size
edo2 <- gseDO(geneList)
# DO (Disease Ontology) Gene Set Enrichment Analysis

dotplot(edo, showCategory=10) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=10) + ggtitle("dotplot for GSEA")

# cnetplot() depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) as a network. GSEA result is also supported with only core enriched genes displayed.
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
p1
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p2
# circular network
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE) 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))

# heatmap
p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2])

# treeplot: hierarchical cluster
edox2 <- pairwise_termsim(edox)
p1 <- treeplot(edox2)
p2 <- treeplot(edox2, hclust_method = "average")
aplot::plot_list(p1, p2, tag_levels='A')

# enrichment map
# mutually overlapping gene sets are tend to cluster together, making it easy to identify functional module
edo <- pairwise_termsim(edo)
p1 <- emapplot(edo)
p2 <- emapplot(edo, cex_category=1.5) # cex_category resize nodes
emapplot(edo,showCategory = 20)
p4 <- emapplot(edo, cex_category=1.5,layout="kk") 
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# biological theme comparison

data(gcSample)
xx <- compareCluster(gcSample, fun="enrichKEGG",
                     organism="hsa", pvalueCutoff=0.05)
# compare multiple clusters of genes
xx <- pairwise_termsim(xx)                     
p1 <- emapplot(xx)
p2 <- emapplot(xx, legend_n=2)
# the number of circles in the bottom left corner can be adjusted

p3 <- emapplot(xx, pie="count")
# proportion of clusters in the pie chart can be adjusted using the pie parameter

p4 <- emapplot(xx, pie="count", cex_category=1.5, layout="kk")
cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# alternative to cnetplot for visualizing the complex association between genes and gene sets. 
# It emphasizes the gene overlapping among different gene sets
upsetplot(edo)

# For GSEA result, it will plot the fold change distributions of different categories 
# (e.g. unique to pathway, overlaps among different pathways)
upsetplot(kk2)

# ridgeplot will visualize expression distributions of core enriched genes for GSEA enriched categories. 
# It helps users to interpret up/down-regulated pathways
ridgeplot(edo2)

# running score and preranked list of GSEA result
gseaplot(edo2, geneSetID = 1, by = "runningScore", title = edo2$Description[1])
gseaplot(edo2, geneSetID = 1, by = "preranked", title = edo2$Description[1])
gseaplot(edo2, geneSetID = 1, title = edo2$Description[1])

# another way to plot GSEA
gseaplot2(edo2, geneSetID = 1, title = edo2$Description[1])

# plot multiple gene sets
gseaplot2(edo2, geneSetID = 1:3)

# display p value table
gseaplot2(edo2, geneSetID = 1:3, pvalue_table = TRUE,
          color = c("#E495A5", "#86B875", "#7DB0DD"), ES_geom = "dot")

# pubmed trend of enriched items
# install.packages("europepmc")
library(cowplot)
terms <- edo$Description[1:5]
p <- pmcplot(terms, 2010:2020) # take a long time...
p2 <- pmcplot(terms, 2010:2020, proportion=FALSE)
plot_grid(p, p2, ncol=2)
