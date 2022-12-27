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

# read and tidy raw data ----------

countsData <- read_csv('../data/rna_seq_counts.csv', col_types = c(entrez = 'c'))

countMatrix <- as.matrix(countsData[2:4])
rownames(countMatrix) <- countsData$entrez

col <- data.frame(sample = c("GFP", "Cap1", "Cap2"),
                  labels = c("control", "sample", "sample"))

# create DESeq dataset -----------
dds <- DESeqDataSetFromMatrix(
  countData = countMatrix,
  colData = col,
  design = ~ labels, # put condition of interests at the end of design formula
  tidy = FALSE
)

# pre-filtering data of too low counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep, ]

# relevel to specify the control group (reference)
dds$labels <- relevel(dds$labels, ref = "control")

# write filtered counts data to disk
counts <- counts(dds)
write.csv(counts, "../data/DESeq2Counts.csv")

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
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
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
  main = "Poisson distances"
)

# plot PCA
plotPCA(vsd, intgroup = c("labels", "sample"))

# differential analysis and result table------
dds <- DESeq(dds)
res <- results(dds)
summary(res) # print a summary of DESeq results


# order our results table by the smallest p value
resOrdered <- res[order(res$padj), ]

# transform DESeqResults object to dataframe for output
resdata <- merge(
  as.data.frame(resOrdered),
  as.data.frame(counts(dds,
                       normalized = TRUE)),
  by = "row.names",
  sort = FALSE)

names(resdata)[1] <- "entrez"

resdata <- left_join(resdata, countsData[c(1, 5)])

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
write_csv(resdata, "../results/DESeq2_DEG.csv")

# volcano plot -------------
resdata %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = type)) +
  geom_hline(aes(yintercept = 2), linetype = "dashed") +
  geom_vline(aes(xintercept = -1), linetype = "dashed") +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  scale_color_discrete(labels = c('Downregulate by Capsid', 'NS', 'Upregulate by Capsid'),
                       type = c('blue', 'grey', 'red'))

fppDEG <- resdata %>%
  filter(type != 0)

# heatmap
fppMatrix <- as.matrix(fppDEG[8:10])

rownames(fppMatrix) <- fppDEG$symbol

pheatmap(fppMatrix,
         scale = 'row',
         #cluster_cols = FALSE,
         angle_col = 0)

# heatmap of IFN-g pathway
ifnDESeq <- resdata %>%
  filter(symbol %in% ifngPathway)

ifnMatrix <- as.matrix(ifnDESeq[8:10])

rownames(ifnMatrix) <- ifnDESeq$symbol

pheatmap(ifnMatrix,
         scale = 'row',
         #cluster_cols = FALSE,
         angle_col = 0)
