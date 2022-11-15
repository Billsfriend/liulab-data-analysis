# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
# official vignette of DESeq2
# follow workflow after salmon and tximeta
# BiocManager::install('vsn')
# BiocManager::install('apeglm')

library(DESeq2)
library(vsn)
library(pheatmap)
library(PoiClaClu) # calculating Poisson distance between samples
library(AnnotationDbi)
library(org.Hs.eg.db)
library(apeglm)
library(tidyverse)
library(ggprism)

# read and tidy raw data ----------
pair1 <- read_csv("data/pair1.csv") %>%
  select(Accession, AbundancesSample, AbundancesControl)

pair2 <- read_csv("data/pair2.csv") %>%
  filter(!is.na(Checked)) %>%
  select(Accession, AbundancesSample, AbundancesControl)

pair3 <- read_csv("data/pair3.csv") %>%
  filter(!is.na(Checked)) %>%
  select(Accession, AbundancesSample, AbundancesControl)

colnames(pair1) <- c("Accession", "BD1", "HC1")
colnames(pair2) <- c("Accession", "BD2", "HC2")
colnames(pair3) <- c("Accession", "BD3", "HC3")

# include all protein detected
allpair <- full_join(pair1, pair2) %>%
  full_join(pair3)

# find protein presenting in all 3 pairs
allAccess <- semi_join(pair1, pair2, by = "Accession") %>%
  semi_join(y = pair3, by = "Accession") %>%
  arrange(Accession)

filtered2 <- filter(pair2, Accession %in% allAccess$Accession) %>%
  arrange(Accession)

filtered3 <- filter(pair3, Accession %in% allAccess$Accession) %>%
  arrange(Accession)

allpair <- left_join(filtered3, filtered2, by = c("Accession", "Description"))

allpair <- left_join(allpair, allAccess, by = c("Accession", "Description"))

write_csv(allpair, "data/bd-neutro-ms.csv")

allpair <- read_csv("data/bd-neutro-ms.csv")

# fill NA with 0
allpair %>%
  replace_na(list(
    BD1 = 0,
    BD2 = 0,
    BD3 = 0,
    HC1 = 0,
    HC2 = 0,
    HC3 = 0
  )) ->
allpair

# abundance * 10 = pseudocount
pseudoCounts <- allpair[-1] * 10
allpair <- bind_cols(allpair[1], pseudoCounts)

count <- as.matrix(dplyr::select(allpair, c(BD1, BD2, BD3, HC1, HC2, HC3)))

rownames(count) <- allpair$Accession

col <- data.frame(sample = c("BD1", "BD2", "BD3", "HC1", "HC2", "HC3"),
                  labels = c("sample", "sample", "sample", "control", "control", "control"),
                  batch = c('1','2','3','1','2','3'))

# create DESeq dataset -----------
dds <- DESeqDataSetFromMatrix(
  countData = count,
  colData = col,
  design = ~ batch + labels, # put condition of interests at the end of design formula
  tidy = FALSE
)

# pre-filtering data of too low counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]

# relevel to specify the control group (reference)
dds$labels <- relevel(dds$labels, ref = "control")

# write filtered counts data to disk
counts <- counts(dds)
write.csv(counts, "data/DESeq2Counts.csv")

# technical replicates (the same biosample) can be collapsed by collapseReplicates command
# biological replicates should not

# transformation below are aim to compress effect of low count genes on distance calculation and PCA plot
# regularized-logarithm transformation. transform counts data to log2
rld <- rlog(object = dds, blind = F)
# or use variance stabilizing transformation instead of rlog:
vsd <- vst(object = dds, blind = 1)
# to perform 'unsupervised transformation' set blind to 1 (default)

# calculate euclidean distance between samples
sampleDists <- dist(t(assay(rld)))

# heatplot of Euclidean distances between samples
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$sample, vsd$labels, sep = " - ")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors,
  main = "Euclidean distances"
)

# use raw counts to calculate Poisson distance between samples
poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix(poisd$dd)
rownames(samplePoisDistMatrix) <- paste(dds$sample, dds$labels, sep = " - ")
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
  clustering_distance_rows = poisd$dd,
  clustering_distance_cols = poisd$dd,
  col = colors,
  main = "Poisson distances"
)

# plot PCA
plotPCA(vsd, intgroup = c("labels", "sample"))

# differential analysis and result table------
dds <- DESeq(dds)
res <- results(dds)
summary(res) # print a summary of DESeq results

# map ENSEMBL id to gene symbol
ens.str <- substr(rownames(res), 1, 18)
res$symbol <- mapIds(org.Hs.eg.db,
  keys = ens.str,
  column = "SYMBOL",
  keytype = "UNIPROT",
  multiVals = "first"
)
res$entrez <- mapIds(org.Hs.eg.db,
  keys = ens.str,
  column = "ENTREZID",
  keytype = "UNIPROT",
  multiVals = "first"
)


# order our results table by the smallest p value
resOrdered <- res[order(res$padj), ]

# transform DESeqResults object to dataframe for output
resdata <- merge(
  as.data.frame(resOrdered),
  as.data.frame(counts(dds,
                       normalized = TRUE)),
  by = "row.names",
  sort = FALSE)

names(resdata)[1] <- "Uniprot"

# visualization ---------
# MA plot
plotMA(res, ylim = c(-5, 5))

# highlight SDEG that p < 0.01 and log2FC > 1
resdata <- resdata %>%
  filter(!is.na(pvalue)) %>%
  mutate(type = as.character(
    (pvalue < 0.01) * (abs(log2FoldChange) > 1)
    * log2FoldChange / abs(log2FoldChange)))

# write DEG data to disk
write.table(resdata,
            file = "results/DEG_all.csv",
            sep = ",",
            row.names = F)

# volcano plot -------------
resdata %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = type)) +
  geom_hline(aes(yintercept = 2), linetype = "dashed") +
  geom_vline(aes(xintercept = -1), linetype = "dashed") +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  scale_color_discrete(labels = c('Enrich in HC', 'NS', 'Enrich in BD'),
                       type = c('blue', 'grey', 'red'))

# annotate unmap uniprot entry
unmap <- read_delim('data/uniprot-accession.txt')
colnames(unmap) <- c('Uniprot', 'symbol')

unmapDESeq <- filter(resdata, is.na(symbol)) %>%
  dplyr::select(-symbol)

unmapDESeq <- left_join(unmapDESeq, unmap)

resdata %>%
  filter(!is.na(symbol)) %>%
  bind_rows(unmapDESeq) ->
  resdata

# genes of interest in FPP metabolism ------
fppList <- read_csv("data/fpp_genes.txt", col_names = "symbol")

fppList$fpp <- TRUE

resdata <- left_join(resdata, fppList)

fppDEG <- filter(resdata, fpp == TRUE)

# volcano plot
resdata %>%
  mutate(fppsize = 1.05 - is.na(fpp)) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = fpp, size = fppsize, alpha = fppsize), position = 'jitter') +
  geom_hline(aes(yintercept = 2), linetype = "dashed") +
  geom_vline(aes(xintercept = -1), linetype = "dashed") +
  geom_vline(aes(xintercept = 1), linetype = "dashed")

# heatmap
fppMatrix <- as.matrix(fppDEG[10:15])

rownames(fppMatrix) <- fppDEG$symbol

pheatmap(fppMatrix[1:6, 1:6],
         scale = 'row',
         cluster_cols = FALSE,
         angle_col = 0)

# barplot of log2FC and pvalue
fppDEG[1:6, 1:17] %>%
  ggplot(aes(x = symbol, y = log2FoldChange)) +
  geom_col(aes(fill = pvalue)) +
  scale_fill_viridis_c(direction = -1) +
  theme_classic(base_size = 16)

# filter out significant DEGs
prettyData <- resdata %>%
  filter(type != "0" & !is.na(symbol)) %>%
  dplyr::select(symbol, HC1, HC2, HC3, BD1, BD2, BD3, type) %>%
  arrange(type)

# plot heatmap of all SDEG
prettyMat <- as.matrix(prettyData[2:7])

rownames(prettyMat) <- prettyData$symbol

pheatmap(prettyMat,
         angle_col = 0,
         scale = 'row',
         show_rownames = FALSE)

# other cytokines in interest -----
cytokine_list <- c('TNFA', 'IL6', 'IL8', 'IL1B', 'IL18', 'IL12A')

resdata %>%
  filter(symbol %in% cytokine_list) -> t

pivot_longer(t, cols = 10:15, names_to = 'sample', values_to = 'normalized_expr') -> t

t$fpp <- c(rep('BD', 3), rep('HC', 3))

df_p_val <- data.frame(
  group1 = "BD",
  group2 = "HC",
  label = 0.1608,
  y.position = 1200
)

t$group <- t$fpp

ggplot(t, aes(group, normalized_expr))+
  stat_summary(fun = 'mean', geom = 'col', aes(fill = group))+
  geom_point(aes(shape = group))+
  add_pvalue(df_p_val) +
  labs(x = 'group',
       title = 'IL18')
  
# cluster of GO pathway in heatmap
bdOverlap <- read_tsv('results/MS-BD-overlap.tsv')
mean