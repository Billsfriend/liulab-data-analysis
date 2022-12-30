# 2022.3.24
# edit in 2022.9.29

# Load libraries
library(Seurat)
library(harmony)
library(tidyverse)
library(SingleCellExperiment)
library(HDF5Array)
library(DelayedArray)

# load seurat objects
sobj <- readRDS('CRC-I/data/Zhang-Yu-2020/zy2020_tumor10x.rds')

sg_sobj <- read_rds('CRC-I/data/seekgene/crc0926.rds')

# sobj data need not normalize by counts per million but only log
sobj[["RNA"]] <- CreateAssayObject(counts = log(sobj[["RNA"]]@counts + 1))
sobj@misc <- list('normalized_tpm')

sobj@assays$RNA@counts[1:5, 1:5]
sobj@assays$RNA@data[1:5, 1:5]

Idents(sobj) <- 'genotype'

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

write_rds(scM_sobj, 'CRC-I/data/Zhang-Yu-2020/scm_sobj.rds')

read_rds('CRC-I/data/Zhang-Yu-2020/scm_sobj.rds') -> scM_sobj

scM_sobj

scM_sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() ->
  scM_sobj

# process merged seurat-----------
scM_sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() ->
  scM_sobj

DimPlot(scM_sobj, group.by = 'batch')

scM_sobj %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors() %>%
  FindClusters() ->
  scM_sobj

DimPlot(scM_sobj)
DimPlot(scM_sobj, group.by = 'Global_Cluster')
DimPlot(scM_sobj, group.by = 'seurat_clusters', split.by = 'batch')

scM_sobj@meta.data %>%
  ggplot(aes(seurat_clusters, fill = batch))+
  geom_bar(position = 'fill')

sobj$Global_Cluster

scM_sobj$seurat_clusters %>% unique()

# process separately -----------
sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() ->
  sobj

sg_sobj %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() ->
  sg_sobj

# identify cell type by public ref ------------
hpca <- read_rds('CRC-I/ref/HumanPrimaryCellAtlas.rds')

annotate_cell_type_by_singler <- function(srt, ref = hpca, ref_label = hpca$label.main, new_label = 'singler_labels'){
  srt %>%
  as.SingleCellExperiment() %>%
  singleR::SingleR(
    ref = ref,
    labels = ref_label,
    clusters = srt$seurat_clusters) ->
  singler_sce

tibble(seurat_clusters = singler_sce@rownames, singler_labels = singler_sce$pruned.labels) -> tmp

colnames(tmp) <- c('seurat_clusters', new_label)

tmp %>%
  right_join(rownames_to_column(srt@meta.data, 'rownames')) %>%
  column_to_rownames('rownames') ->
  srt@meta.data
}

annotate_cell_type_by_singler(sobj) -> sobj@meta.data

annotate_cell_type_by_singler(sg_sobj) -> sg_sobj@meta.data

annotate_cell_type_by_singler(scM_sobj) -> scM_sobj@meta.data

# identify cell type by 10x Zhang-Yu 2020 as ref ----------
sobj %>% as.SingleCellExperiment() -> sce10x

annotate_cell_type_by_singler(
  sg_sobj,
  ref = sce10x,
  ref_label = sce10x$Sub_Cluster, 
  new_label = 'Sub_Cluster') ->
sg_sobj@meta.data


# save integrated file ---------
write_rds(sg_sobj, 'CRC-I/data/crc0926_singler.rds')
write_rds(sobj, "CRC-I/data/Zhang-Yu-2020/tenx_singler.rds")
write_rds(scM_sobj, 'CRC-I/data/scMerge_CRC221230.rds')