library(Seurat)
library(tidyverse)
library(SingleR)

Read10X('IgG-IgM-B/data/Liu-Qi_2022_LLPC/d14-bm1-alevin/') -> d14_bm1

Read10X_h5('IgG-IgM-B/data/Liu-Qi_2022_LLPC/d14-bm1.h5') -> d14_bm1

CreateSeuratObject(d14_bm1, min.cells = 1, min.features = 1) -> sobj

sobj %>% subset(nCount_RNA >= 200)

VlnPlot(sobj, 'nFeature_RNA')

# alevin output ------------
data.table::fread('IgG-IgM-B/data/Liu-Qi_2022_LLPC/D14-BM1-GE-alevin.gz') -> mtx

column_to_rownames(mtx, 'V1') %>%
  t() %>%
  CreateSeuratObject(min.cells = 1, min.features = 200) ->
  sobj

DropletUtils::write10xCounts(
  'IgG-IgM-B/data/Liu-Qi_2022_LLPC/d14-bm1-alevin/',
  sobj@assays$RNA@counts,
  overwrite = TRUE,
  version = '3')

sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20) ->
  sobj

DimPlot(sobj)

FeaturePlot(sobj, c('Ighg1', 'Cd19'))

# identify per cluster ------------
celldex::ImmGenData() -> singler_ref

per_clust_idented <- sobj %>%
  as.SingleCellExperiment() %>%
  SingleR(
    ref = singler_ref,
    labels = singler_ref$label.fine,
    clusters = sobj$seurat_clusters,
    de.method = "wilcox"
  )

plotScoreHeatmap(per_clust_idented,
                 show.pruned = TRUE,
                 show_colnames = TRUE)

# cell ranger aggr output -----------
Read10X_h5('IgG-IgM-B/data/Liu-Qi_2022_LLPC/LLPC_aggr.h5') -> aggr_mtx

aggr_mtx[1:5,1:5]

read_csv('IgG-IgM-B/data/Liu-Qi_2022_LLPC/LLPC_aggr_table.csv') %>%
  add_column(orig.ident = factor(1:18)) %>%
  select(-molecule_h5) -> aggr_table

tibble(barcodes = colnames(aggr_mtx)) %>%
  separate(1, into = c('cell', 'orig.ident')) %>%
  left_join(aggr_table) %>%
  unite(col = 'new_barcodes', c(cell, sample_id), sep = '-') %>%
  pull('new_barcodes') ->
  colnames(aggr_mtx)

CreateSeuratObject(aggr_mtx,
                   min.cells = 1,
                   min.features = 1,
                   names.delim = '-',
                   names.field = 2) -> sobj

VlnPlot(sobj, 'nCount_RNA')

write_rds(sobj, 'IgG-IgM-B/data/Liu-Qi_2022_LLPC/aggr.rds', 'gz') %>% system.time()

system.time(read_rds('IgG-IgM-B/data/Liu-Qi_2022_LLPC/aggr.rds') -> sobj)