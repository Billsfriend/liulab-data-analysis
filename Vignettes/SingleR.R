# http://bioconductor.org/books/devel/SingleRBook/introduction.html
# BiocManager::install('celldex')
# BiocManager::install('SingleR')
# BiocManager::install('scRNAseq')
library(tidyverse)

library(celldex)
hpca.se <- celldex::HumanPrimaryCellAtlasData() # Human Primary Cell Atlas data as a SummarizedExperiment object containing a matrix of log-expression values with sample-level labels
# these datasets will be used as reference
write_rds(hpca.se, "ref/HumanPrimaryCellAtlas.rds")

library(scRNAseq)
hESCs <- LaMannoBrainData("human-es") # some human embryonic stem cells as example test data
hESCs <- hESCs[, 1:100] # sample 100 cells for quick demo

# use our hpca.se reference to annotate each cell in hESCs via the SingleR() function

library(SingleR)
pred.hesc <- SingleR(
  test = hESCs,
  ref = hpca.se,
  assay.type.test = 1,
  labels = hpca.se$label.main
)
# Labels are shown before fine-tuning (first.labels), after fine-tuning (labels) and after pruning (pruned.labels), along with the associated scores

# Summarizing the distribution:
table(pred.hesc$labels)

# The above example uses SummarizedExperiment objects, but the same functions will accept any (log-)normalized expression matrix

# Using single-cell reference ----------
sceM <- MuraroPancreasData()

# One should normally do cell-based quality control at this point, but for brevity's sake, we will just remove the unlabelled libraries here.
sceM <- sceM[, !is.na(sceM$label)]

# SingleR() expects reference datasets to be normalized and log-transformed.
# BiocManager::install('scuttle')
library(scuttle)
sceM <- logNormCounts(sceM)

# single-cell test data
sceG <- GrunPancreasData()
sceG <- sceG[, colSums(counts(sceG)) > 0] # Remove libraries with no counts.
sceG <- logNormCounts(sceG)
sceG <- sceG[, 1:100] # sample 100 cells

# BiocManager::install('scran')
pred.grun <- SingleR(
  test = sceG,
  ref = sceM,
  labels = sceM$label,
  de.method = "wilcox"
)
# Here, we will use the Wilcoxon ranked sum test to identify the top markers for each pairwise comparison between labels. This is slower but more appropriate for single-cell data compared to the default marker detection algorithm (which may fail for low-coverage data where the median is frequently zero)
table(pred.grun$labels)

# Annotation diagnostics
plotScoreHeatmap(pred.grun)

# Low deltas indicate that the assignment is uncertain, which is especially relevant if the cell's true label does not exist in the reference
plotDeltaDistribution(pred.hesc, ncol = 3)

# pruneScores() function will remove potentially poor-quality or ambiguous assignments based on the deltas
summary(is.na(pred.hesc$pruned.labels))

# another simple diagnostic method is to examine the cluster marker expression of each labelled cells
all.markers <- metadata(pred.hesc)$de.genes
natSE$labels <- pred.hesc$labels

# examine labeling of a specific cell type ----
# BiocManager::install('scater')
library(scater)
# Beta cell-related markers
plotHeatmap(natSE,
  order_columns_by = "labels",
  features = unique(unlist(all.markers$Eosinophils))
)

# other references -----------
# immune references from celldex
immgen <- ImmGenData() # most exhaustive coverage of dizzying number of subtypes, which is all murine

diceRef <- celldex::DatabaseImmuneCellExpressionData() # particularly useful to those interested in CD4+ T cell subsets, though the lack of CD4+ central memory, effector memory, DC and B cell subset sample

NHref <- NovershternHematopoieticData() # likely the best option for bone marrow samples

MIref <- MonacoImmuneData() # human immune reference that best covers all of the bases for a typical PBMC sample

# try on my datasets -----------
library(Seurat)
natTpm <- read_rds("data/nature/CountNat.rds")
head(natTpm@assays$RNA@data)

# transform seurat object into SummarizedExperiment
natSE <- as.SingleCellExperiment(natTpm)

# or just get log-count data matrix
fetchNat <- FetchData(natTpm, vars = "ident", slot = "data")

DataNat <- GetAssayData(natTpm, "data")

library(BiocParallel)
pred.dice <- SingleR(
  test = natSE,
  ref = diceRef,
  de.method = "wilcox",
  # clusters = natSE$seurat_clusters,
  labels = diceRef$label.fine,
  BPPARAM = MulticoreParam()
)

table(pred.hesc$labels)

plotScoreHeatmap(pred.hesc)

pred.hesc <- SingleR(
  test = natSE,
  ref = hpca.se,
  de.method = "wilcox",
  labels = hpca.se$label.fine,
  BPPARAM = MulticoreParam()
)

table(pred.hesc$labels)

# cluster-based annotation
library(scran)
dec <- modelGeneVarByPoisson(natSE)
natSE <- denoisePCA(natSE,
                    dec,
                    subset.row = getTopHVGs(dec, n = 5000))

library(bluster)
colLabels(natSE) <- clusterRows(reducedDim(natSE), NNGraphParam())

library(scater)
set.seed(117)
natSE <- runTSNE(natSE, dimred = "PCA")
plotTSNE(natSE, colour_by = "label")

SingleR(natSE,
        diceRef,
        clusters = colLabels(natSE),
        labels = diceRef$label.main)

diceCluster <- SingleR(natSE,
                       diceRef,
                       clusters = natSE$seurat_clusters,
                       labels = diceRef$label.fine)

table(diceCluster$labels)

plotScoreHeatmap(diceCluster)

plotDeltaDistribution(diceCluster)

MICluster <- SingleR(natSE,
                     MIref,
                     clusters = natSE$seurat_clusters,
                     labels = MIref$label.fine)

table(MICluster$labels)

plotScoreHeatmap(MICluster)

hpcaCluster <- SingleR(natSE,
                       hpca.se,
                       clusters = natSE$seurat_clusters,
                       labels = hpca.se$label.fine)

table(hpcaCluster$labels)

plotScoreHeatmap(hpcaCluster)

# SingleR results can be easily add to seurat object as metadata
natTpm[["singler.labels"]] <- singler.results$labels

# Or if `method="cluster"` was used:
natTpm[["SingleR.MI.labels"]] <-
  MICluster$labels[match(natTpm[[]][["seurat_clusters"]], rownames(MICluster))]

natTpm[["SingleR.DICE.labels"]] <-
  diceCluster$labels[match(natTpm[[]][["seurat_clusters"]], rownames(diceCluster))]

natTpm[["SingleR.MI.labels"]] <-
  MICluster$labels[match(natTpm[[]][["seurat_clusters"]], rownames(MICluster))]

natTpm[["SingleR.hpca.labels"]] <-
  hpcaCluster$labels[match(natTpm[[]][["seurat_clusters"]], rownames(hpcaCluster))]

write_rds(natTpm, "data/nature/CountNat.rds")
