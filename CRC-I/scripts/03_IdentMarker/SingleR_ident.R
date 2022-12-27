library(Seurat)
library(celldex)
library(SingleR)
library(BiocParallel)
library(tidyverse)

ref_huPrCellAtlas <- read_rds("ref/HumanPrimaryCellAtlas.rds")

# read your seurat data here ---------
to_ident <- read_rds('data/Wang2022/harmonySobj.rds')

rownames(to_ident@meta.data) -> to_ident$barcode

# identify per cell on main types -------
per_cell_idented <- to_ident %>%
  as.SingleCellExperiment() %>%
  SingleR(
  ref = ref_huPrCellAtlas,
  labels = ref_huPrCellAtlas$label.main,
  de.method = "wilcox",
  BPPARAM = MulticoreParam()
)

table(per_cell_idented$pruned.labels)

plotScoreHeatmap(per_cell_idented)

summary(is.na(per_cell_idented$pruned.labels))

# identify per cluster ------------
per_clust_idented <- to_ident %>%
  as.SingleCellExperiment() %>%
  SingleR(
    ref = ref_huPrCellAtlas,
    labels = ref_huPrCellAtlas$label.main,
    clusters = to_ident$seurat_clusters,
    de.method = "wilcox",
    BPPARAM = MulticoreParam()
  )

table(per_clust_idented$pruned.labels)
# note: CMP - Common Myeloid Progenitor

plotScoreHeatmap(per_clust_idented,
                 show.pruned = TRUE,
                 show_colnames = TRUE)

summary(is.na(per_clust_idented$pruned.labels))

to_ident[["SingleR.labels"]] <-
  per_clust_idented$labels[match(to_ident[[]][["seurat_clusters"]], rownames(per_clust_idented))]

DimPlot(to_ident,
        group.by = 'SingleR.labels',
        split.by = 'genotype') +
  scale_color_brewer(type = 'qual', palette = 'Set1') +
  labs(title = 'UMAP of cell types in CRC')

write_rds(to_ident, 'data/Wang2022/singleR_main.rds')

# identify immune cells per cluster
to_immu_ident <- subset(to_ident, SingleR.labels %in% c('B_cell', 'T_cells', 'Macrophage'))

DimPlot(to_immu_ident)

write_rds(to_immu_ident, 'data/Wang2022/singleR_immu.rds')
