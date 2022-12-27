# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#countmat
# official vignette of DESeq2
# follow workflow after salmon and tximeta
# BiocManager::install('vsn')
# BiocManager::install('apeglm')
# try to work with snakemake

library(DESeq2)
library(tidyverse)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(apeglm)
library(BiocParallel)
register(MulticoreParam(4)) # add argument parallel=TRUE to time-consuming function to speed up

gse <- read_rds('data/geneSumExp.rds')

# force condition character as factor so that can be used in design
gse$condition <- as.factor(gse$condition)

dds <- DESeqDataSet(gse, design = ~ condition)

# pre-filtering data of too low counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

# relevel to specify the control group (reference)
levels(dds$condition) <- c('IgG', 'IgM')
dds$condition <- relevel(dds$condition, ref = "IgM")

# write filtered counts data to disk
counts <- counts(dds)
counts <- as_tibble(counts)
counts$symbol <- dds@rowRanges$gene_name
counts$ensembl <- dds@rowRanges$gene_id
write.csv(counts, 'data/Counts.csv')

# technical replicates (the same biosample) can be collapsed by collapseReplicates command
# biological replicates should not

# transformation below are aim to compress effect of low count genes on distance calculation and PCA plot
# regularized-logarithm transformation. transform counts data to log2
# rld <- rlog(object=dds,blind=F)
# rlog need replicates, which is not applicable in this sample

# or use variance stabilizing transformation instead of rlog:
vsd <- vst(object=dds,blind=1)
# to perform 'unsupervised transformation' set blind to 1 (default) 

dds <- DESeq(dds)
res <- results(dds)
summary(res) # print a summary of DESeq results

# map ENSEMBL id to gene symbol
ens.str <- substr(rownames(res), 1, 18)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
ressym <- subset(res, !is.na(res$symbol)) # some ensembl id have no corresponding symbol...
rownames(ressym) <- ressym$symbol

# order our results table by the smallest p value
resOrdered <- res[order(res$padj),]

# lower the FDR threshold from default 10% to 5%
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)

# set a log2FC threshold of 1
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

# transform DESeqResults object to dataframe for output
resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata)[1] <- "Gene"

# select only genes of p.adj < 0.05 and log2FC > 1
resdata_filter <- subset(resdata,resdata$padj<0.05 & abs(resdata$log2FoldChange) >1)

# write DEG data to disk
write.table(resdata_filter, file="DiffAnalysis/DEG_SalTx_0.05.csv",sep=",",row.names=F)

write.table(resdata,file="DiffAnalysis/DEG_SalTx_results.all.csv",sep=",",row.names=F)

# plot most differentiated gene expression
topGene <- rownames(res)[which.min(res$padj)]
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("condition","names"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = condition, y = count, color = names)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

# Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. e.g. MAplot
resLFC <- lfcShrink(dds, coef = 'condition_cko_vs_con',type="apeglm")
plotMA(resLFC, ylim = c(-5, 5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
}) # highlight a top gene in MA plot

# plot most variable genes in heatmap
library("genefilter")
topVarGenes <- head(order(rowVars(assay(ressym)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
# show relative VST-transformed values across samples
anno <- as.data.frame(colData(vsd)[, c("names","condition")])
pheatmap(mat, annotation_col = anno)
