library(tidyverse)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(pheatmap)

# reading in the DEG list ----------
resdata <- read_csv('results/DEG_all.csv') %>%
  filter(!is.na(symbol))

# a list of all genes in DEG analysis
unilist <- resdata$symbol

uniEntrez <- resdata$entrez

# a numeric list of log2fc for GSEA
foldchange <- resdata$log2FoldChange
names(foldchange) <- unilist

foldchange <- sort(foldchange, decreasing = TRUE)

BDdeg <- filter(resdata, type == '1')
HCdeg <- filter(resdata, type == '-1')

# looser criteria for DEGs
BDdeg <- subset(resdata, pvalue < 0.05 & log2FoldChange > 1)

# pathway enrichment analysis --------

# see overlaps in GO:BP sets and make no significance cutoff
ggoHC <- groupGO(gene = HCdeg$symbol,
                 keyType = 'SYMBOL',
                 OrgDb = org.Hs.eg.db,
                 ont = 'BP',
                 level = 8)

t <- as_tibble(ggoHC@result)

# over-represent analysis on GO
oraGoHC <- enrichGO(gene = HCdeg$symbol,
                  keyType = 'SYMBOL',
                  OrgDb = org.Hs.eg.db,
                  universe = unilist,
                  ont = 'ALL',
                  minGSSize = 3,
                  pool = TRUE)

head(oraGoHC)

oraGoBD <- enrichGO(gene = BDdeg$symbol,
                    keyType = 'SYMBOL',
                    OrgDb = org.Hs.eg.db,
                    universe = unilist,
                    ont = 'BP',
                    minGSSize = 3)

oraKeggBD <- enrichKEGG(gene = BDdeg$entrez,
                        universe = uniEntrez,
                        minGSSize = 3,
                        pvalueCutoff = 0.1,
                        qvalueCutoff = 0.5)

# GSEA
gsaGoBD <- gseGO(geneList = foldchange,
                 OrgDb = org.Hs.eg.db,
                 ont = 'ALL',
                 keyType = 'SYMBOL')

# visualization ----------
dotplot(oraGoHC)

dotplot(oraGoBD)

heatplot(oraGo, showCategory = 10, foldChange = foldchange)
# why no foldchange color?

# extract symbols in enriched pathway
hcEnrich <- oraGo@result[["geneID"]][1]

hcDEG <- str_split(hcEnrich, '/', simplify = T)

resdata %>%
  filter(symbol %in% hcDEG) ->
  hcPathway

hcMatrix <- as.matrix(hcPathway[10:15])

rownames(hcMatrix) <- hcPathway$symbol
colnames(hcMatrix) <- c('BD1','BD2','BD3','HC1','HC2','HC3')

pheatmap(hcMatrix,
         scale = 'row',
         main = 'GO:CC extracellular vesicle',
         angle_col = 0)

# re-draw heatmap of all degs
resdata %>%
  filter(type != 0) ->
  allDeg

degMatrix <- as.matrix(allDeg[10:15])

rownames(degMatrix) <- allDeg$symbol
colnames(degMatrix) <- c('BD1','BD2','BD3','HC1','HC2','HC3')

# add attribute to some cases by filter-add_column-right_join
pathwayAnnotation <- allDeg %>%
  dplyr::filter(symbol %in% hcDEG) %>%
  add_column(pathway = 'Extracellular vesicle') %>%
  right_join(allDeg)

# row annotation data frame
matrixAnno <- as.data.frame(pathwayAnnotation[17])
rownames(matrixAnno) <- pathwayAnnotation$symbol

# whole heatmap
pheatmap(degMatrix,
         scale = 'row',
         main = 'All DEGs in HC-BD neutrophils',
         angle_col = 0,
         annotation_row = matrixAnno,
         annotation_legend = FALSE)

# msigdb overlap
overlapBD <- read_tsv('results/MS-BD-overlap.tsv')

overlapBD %>%
  ggplot(aes(x = Gene_Set, y = Counts)) +
  geom_col(aes(fill = `FDR-value`)) +
  coord_flip() +
  scale_fill_viridis_c(direction = -1)

overlapHC <- read_tsv('results/MS-HC-overlap.tsv')

overlapHC %>%
  ggplot(aes(x = Gene_Set, y = Counts)) +
  geom_col(aes(fill = `FDR-value`)) +
  coord_flip() +
  scale_fill_viridis_c(direction = -1)

# heatmap of pathway
extracell_pathway <- read_csv('results/extra_vesicle_pathway.txt', col_names = FALSE)

resdata %>%
  filter(symbol %in% extracell_pathway$X1) %>%
  dplyr::select(c(HC1, HC2, HC3, BD1, BD2, BD3)) %>%
  pheatmap(scale = 'row')

hc_pathway <- read_tsv('results/MS-HC-pathwayGene.txt')

hc_pathway %>%
  filter(!is.na(GOCC_SECRETORY_VESICLE)) %>%
  pull(`Gene Symbol`)
