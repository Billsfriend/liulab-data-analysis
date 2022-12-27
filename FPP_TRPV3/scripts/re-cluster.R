library(Seurat)
library(harmony)
library(tidyverse)
library(celldex)
library(SingleR)
library(BiocParallel)
library(Platypus)

sobj <- read_rds('data/Wang2022/singleR_immu.rds')

metadata <- sobj@meta.data

metadata$genotype <- NA

metadata$genotype[which(metadata$patient == 'P5')] <- 'TT' 

metadata$genotype[which(metadata$patient %in% c('P3','P4','P7','P9'))] <- 'IT' 

metadata$genotype[which(metadata$patient %in% c('P1','P2','P6','P8'))] <- 'II' 

metadata <- unite(metadata,
                  patient,
                  tissue,
                  type,
                  col = 'batch',
                  sep = '_',
                  remove = FALSE)

sobj@meta.data <- metadata

# re-clustering ---------
sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() ->
  sobj

DimPlot(sobj,
        reduction = 'pca',
        group.by = 'patient')

harmonySobj <- RunHarmony(sobj,
                          'batch',
                          max.iter.harmony = 20)

DimPlot(harmonySobj,
        reduction = 'harmony',
        group.by = 'patient')

# run UMAP and clustering
harmonySobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity() ->
  harmonySobj

DimPlot(harmonySobj, label = TRUE)

DimPlot(harmonySobj, reduction = 'harmony')

RunTSNE(harmonySobj) -> harmonySobj

DimPlot(harmonySobj, reduction = 'tsne')

# identify cluster
allMarkers <- FindAllMarkers(harmonySobj,
                             only.pos = TRUE)

allMarkers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = p_val_adj) ->
  topMarkers

write_csv(topMarkers, 'results/cluster/wang2022immuneTopMarkers.csv')

CD4TandOtherImmune <- celldex::DatabaseImmuneCellExpressionData()

PbmcImmune <- celldex::MonacoImmuneData()

# identify per cluster ------------
per_clust_idented <- harmonySobj %>%
  as.SingleCellExperiment() %>%
  SingleR(
    ref = CD4TandOtherImmune,
    labels = CD4TandOtherImmune$label.main,
    clusters = harmonySobj$seurat_clusters,
    de.method = "wilcox",
    BPPARAM = MulticoreParam()
  )

table(per_clust_idented$pruned.labels)
# note: CMP - Common Myeloid Progenitor

plotScoreHeatmap(per_clust_idented,
                 show_colnames = TRUE,
                 order.by = 'clusters',
                 clusters = 0:6)

harmonySobj$cell_type <- NA

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 0)] <- 'B_cells'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 1)] <- 'Macrophages'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 2)] <- 'CD4+T_cells'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 3)] <- 'B_cells'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 4)] <- 'DCs'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 5)] <- 'CD8+T_cells'

harmonySobj$cell_type[which(harmonySobj$seurat_clusters == 6)] <- 'Monocytes'

DimPlot(harmonySobj, group.by = 'cell_type')

write_rds(harmonySobj, 'data/Wang2022/identImmune.rds')

write.csv(harmonySobj@meta.data, 'data/Wang2022/identImmune.csv')

# differential analysis ---------
harmonySobj <- read_rds('data/Wang2022/identImmune.rds')

# B cell marker
FeaturePlot(harmonySobj, c('CD79A'))

# plasma cell marker
FeaturePlot(harmonySobj, c('MIF', 'XBP1', 'CD27', 'MZB1'))

# T cell marker
FeaturePlot(harmonySobj, c('CD3E', 'IL7R'))

# split by cell type
typeList <- SplitObject(harmonySobj, split.by = 'cell_type')

findmarkInTT <- function(srt){try(FindMarkers(
  srt,
  group.by = 'genotype',
  ident.1 = 'TT',
  fc.name = "avg_logFC"))
}

degList <- lapply(typeList, findmarkInTT)

degList$Fibroblasts %>%
  GEX_volcano(input.type = 'findmarkers',
              condition.1 = 'TT',
              condition.2 = 'II+IT',
              n.label.up = 6)

typeList <- SplitObject(to_ident, split.by = 'SingleR.labels')
