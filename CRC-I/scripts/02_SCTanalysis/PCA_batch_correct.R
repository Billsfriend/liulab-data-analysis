# 2022.3.24
# edit in 2022.9.29

# Load libraries
library(Seurat)
library(harmony)
library(tidyverse)
library(SingleCellExperiment)

# load seurat objects
sobj <- readRDS('CRC-I/data/Zhang-Yu-2020/zy2020_tumor10x.rds')

sg_sobj <- read_rds('CRC-I/data/seekgene/crc0926.rds')

# sobj data need not normalize by counts per million but only log
sobj[["RNA"]] <- CreateAssayObject(counts = log1p(sobj[["RNA"]]@counts))
sobj@misc <- list('normalized_tpm')

sobj@assays$RNA@counts[1:5, 1:5]
sobj@assays$RNA@data[1:5, 1:5]

sg_sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() -> sg_sobj

sg_sobj@assays$RNA@counts[1:10, 1:10]
sg_sobj@assays$RNA@data[1:10, 1:10]

# visualize CD45 expression as an ref to see if expression scale is the same between 2 data sets.
merge(sobj, sg_sobj) -> mrg_obj

mrg_obj@meta.data  %>%
  mutate(batch = case_when(
    is.na(Platform) ~ 'seekgene',
    TRUE ~ '10x'
  )) -> 
  mrg_obj@meta.data

VlnPlot(mrg_obj, c('PTPRC','ACTB'), group.by = 'batch') &
  xlab('Platform')

# try scMerge ---------
data("segList", package = "scMerge")

system.time(mrg_obj %>%
  as.SingleCellExperiment() %>%
  scMerge::scMerge(ctl = segList$human$human_scSEG,
                   assay_name = 'scMerge',
                   kmeansK = c(21, 21),
                   BSPARAM = BiocSingular::RandomParam(),
                   BPPARAM = BiocParallel::MulticoreParam(workers = 6),
                   verbose = TRUE) -> scm_sce)

scm_sce@assays@data$scMerge[1:5,1:5]

sobj$Global_Cluster

# process separately -----------
sobj@meta.data %>%
  select(contains(c('Cluster', 'genotype', 'sample'))) ->
  sobj@meta.data

sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() ->
  sobj

DimPlot(sobj, group.by = 'seurat_clusters', reduction = 'pca')
DimPlot(sobj, group.by = 'Sub_Cluster')

sg_sobj %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() ->
  sg_sobj

# identify cell type by public ref ------------
hpca <- read_rds('CRC-I/ref/HumanPrimaryCellAtlas.rds')

monaco <- celldex::MonacoImmuneData()

annotate_cell_type_by_singler <- function(
    srt,
    ref = hpca,
    ref_label = hpca$label.main,
    new_label = 'singler_labels',
    de_method = 'classic'){
  srt %>%
  as.SingleCellExperiment() %>%
  SingleR::SingleR(
    ref = ref,
    labels = ref_label,
    clusters = srt$seurat_clusters,
    de.method = de_method) ->
  singler_sce

  SingleR::plotScoreHeatmap(singler_sce, 
                            show.pruned = TRUE,
                            show_colnames = TRUE)
  
  tibble(seurat_clusters = singler_sce@rownames, singler_labels = singler_sce$pruned.labels) %>%
  rename_with(~new_label, singler_labels) %>%
  right_join(rownames_to_column(srt@meta.data, 'rownames')) %>%
  column_to_rownames('rownames')
}

annotate_cell_type_by_singler(sobj) -> sobj@meta.data

DimPlot(sobj)

annotate_cell_type_by_singler(sg_sobj) -> sg_sobj@meta.data

DimPlot(sg_sobj)
DimPlot(sg_sobj, group.by = 'singler_labels')

# check_duplicates is a latent parameter in RunTSNE, set to FALSE if error is thrown
sg_sobj <- RunTSNE(sg_sobj, check_duplicates = FALSE)

TSNEPlot(sg_sobj, group.by = 'singler_labels')

# remove stromal cells and re-analyze crc0926 ----
sg_sobj <- subset(sg_sobj, singler_labels %in% c(
  'B_cell',
  'DC',
  'Monocyte',
  'NK_cell',
  'T_cells'
))

sg_sobj <- sg_sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters()

DimPlot(sg_sobj)

sg_sobj@meta.data <- annotate_cell_type_by_singler(
  sg_sobj,
  ref = monaco,
  ref_label = monaco$label.fine,
  new_label = 'monaco_label')

DimPlot(sg_sobj, group.by = 'monaco_label')

sg_sobj@meta.data %>%
  dplyr::count(seurat_clusters, monaco_label)

tribble(~seurat_clusters, ~eco_label,
       6 , 'B.cells',
       7 , 'Monocytes.and.Macrophages',
       11 , 'Monocytes.and.Macrophages',
       1 , 'PMNs',
       12 , 'Dendritic.cells',
       8 , 'NK.cells',
       2 , 'PCs',
       9 , 'PCs',
       10 , 'PCs',
       13 , 'PCs',
       14 , 'Dendritic.cells',
       0 , 'CD4.T.cells',
       3 , 'CD8.T.cells',
       4 , 'CD8.T.cells',
       5 , 'CD4.T.cells') -> eco_tbl

sg_sobj@meta.data <- sg_sobj@meta.data %>%
  rownames_to_column('rownames') %>%
  mutate(seurat_clusters = as.numeric(seurat_clusters)) %>%
  left_join(eco_tbl) %>%
  column_to_rownames('rownames')

sg_sobj$genotype <- 'TT'

eco_input <- sg_sobj@meta.data %>%
  rownames_to_column('ID') %>%
  select(c(ID, eco_label, genotype))

colnames(eco_input) <- c('ID', 'CellType', 'genotype')

write_tsv(eco_input, 'CRC-I/results/crc0926_immune_eco_anno.tsv')

sg_sobj@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column('Cell-ID') -> df

write_tsv(df, 'CRC-I/results/crc0926_immune_mat.tsv')

sobj@assays$RNA@data %>%
  as.data.frame() %>%
  rownames_to_column('Cell-ID') %>%
  write_tsv('CRC-I/results/zy2020_tme_mat.tsv')

sobj@meta.data %>%
  rownames_to_column('ID') %>%
  mutate(CellType = case_when(
    str_detect(Sub_Cluster, 'CD4') ~ 'CD4.T.cells',
    str_detect(Sub_Cluster, 'CD8') ~ 'CD8.T.cells',
    str_detect(Sub_Cluster, 'NK') ~ 'NK.cells',
    str_detect(Sub_Cluster, 'Plasma') ~ 'PCs',
    str_detect(Sub_Cluster, 'hB') ~ 'B.cells',
    str_detect(Sub_Cluster, 'DC') ~ 'Dendritic.cells',
    str_detect(Sub_Cluster, 'Mast') ~ 'Mast.cells',
    str_detect(Sub_Cluster, 'hM') ~ 'Monocytes.and.Macrophages',
    TRUE ~ Sub_Cluster
  )) %>%
  select(c(ID, CellType, genotype)) %>%
  write_tsv('CRC-I/results/zy2020_immune_eco_anno.tsv')
  
# identify cell type by 10x Zhang-Yu 2020 as ref ----------
sobj %>% as.SingleCellExperiment() -> sce10x

new_meta <- annotate_cell_type_by_singler(
  sg_sobj,
  ref = sce10x,
  ref_label = sce10x$Sub_Cluster, 
  new_label = 'Sub_Cluster',
  de_method = 'wilcox')

# save integrated file ---------
write_rds(sg_sobj, 'CRC-I/data/crc0926_singler.rds')
write_rds(sobj, "CRC-I/data/Zhang-Yu-2020/tenx_singler.rds")
write_rds(scM_sobj, 'CRC-I/data/scMerge_CRC221230.rds')