library(Seurat)
library(tidyverse)

# read in matrix data ----------
data.table::fread('data/GSE160269_UMI_matrix_Fibroblast.txt.gz') -> count_mtx

count_mtx[1:5, 1:5]

# check if rownames are unique
count_mtx <- column_to_rownames(count_mtx, 'V1')

sobj <- CreateSeuratObject(count_mtx,
                           names.delim = '-',
                           min.cells = 1,
                           min.features = 1)

# modify metadata ----------
metadata <- sobj@meta.data

metadata <- rownames_to_column(metadata, 'barcode')

ITlist <- read.table('data/ITlist.txt')

metadata %>%
  mutate(genotype = case_when(
    orig.ident %in% c('P8T', 'P63T', 'P130T') ~ 'TT',
    orig.ident %in% ITlist$V1 ~ 'IT',
    orig.ident %in% c("P126N", "P127N", "P128N", "P130N") ~ 'Normal',
    TRUE ~ 'II'
  )) %>%
  separate(barcode, into = c('foo', 'bar','tag'), remove = 0) -> flt_meta

flt_meta %>%
  mutate(barcode2 = case_when(
    str_detect(foo, '126|127|128|130') ~ barcode,
    TRUE ~ str_remove(barcode, 'T')
  )) -> flt_meta

# read from excel ------------
meta_from_paper <- readxl::excel_sheets('data/source_ESM.xlsx') %>%
  set_names() %>%
  map(read_excel,
      path = 'data/source_ESM.xlsx',
      skip = 1)

pub_cell_meta <- meta_from_paper[["Metadata for fibroblasts"]] %>%
  select(c('cell', 'celltype', 'tissue', 'sample')) %>%
  mutate(barcode2 = str_replace_all(cell, '\\.', '-'))

pub_cell_meta %>%
  separate(barcode, into = c('foo', 'bar','tag'), remove = 0) -> pub_cell_meta

left_join(flt_meta, pub_cell_meta, by = 'barcode2') -> bind_meta

bind_meta %>%
  column_to_rownames('barcode') -> sobj@meta.data

sobj <- subset(sobj, genotype != 'Normal')

sobj <- SetIdent(sobj, value = 'genotype')

write_rds(sobj, 'data/zhang-Lin_fibroblast.rds')

VlnPlot(sobj, 'nFeature_RNA', pt.size = 0)

# too many cells to handle on my PC
# sample 4000 cells for II and IT
metadata %>%
  filter(genotype == 'II') %>%
  sample_n(size = 4000) %>%
  rownames() ->
  IIsample

metadata %>%
  filter(genotype == 'IT') %>%
  sample_n(size = 4000) %>%
  rownames() ->
  ITsample

metadata %>%
  filter(genotype == 'TT') %>%
  rownames() ->
  TTsample

sobj <- subset(sobj,
               cells = c(IIsample, ITsample, TTsample))

write_rds(sobj, 'data/zhang_T_sample.rds')

# use runSrublet to remove doublet ------
find_doublet_seurat <- function(s_obj){
  s_obj %>%
    as.SingleCellExperiment() %>%
    singleCellTK::runScrublet() %>%
    as.Seurat()
}

sobj %>%
  SplitObject(split.by = 'orig.ident') %>%
  map(safely(find_doublet_seurat)) %>%
  flatten() %>%
  compact() %>%
  discard(is.list) %>%
  reduce(merge) %>%
  subset(scrublet_call == 'Singlet') ->
  sobj

# pre-process --------
sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() ->
  sobj



DEenrichRPlot(sobj, 
              ident.1 = 'NAF2',
              enrich.database = 'GO_Biological_Process_2021',
              max.genes = 3000,
              return.gene.list = 1) -> delist

DimPlot(sobj, group.by = 'genotype')

# use harmony to clear batch effect ------
harmonySobj <- RunHarmony(sobj, 
                          'genotype')

DimPlot(harmonySobj,
        reduction = 'harmony')

# run UMAP and clustering
harmonySobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.6) %>% 
  identity() ->
  harmonySobj

DimPlot(harmonySobj, reduction = 'umap')+
  scale_color_viridis_d(option = 'turbo', begin = 0.1)
DimPlot(harmonySobj, 
        reduction = 'umap',
        split.by = 'genotype')

RunTSNE(harmonySobj, reduction = 'harmony') -> harmonySobj
TSNEPlot(harmonySobj) +
  scale_color_viridis_d(option = 'turbo', begin = 0.1)

# find markers for clusters
harmonySobj %>%
  FindAllMarkers(only.pos = TRUE,
                 min.pct = 0.1) %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = avg_log2FC) ->
  topMarkers

write_csv(topMarkers, 'results/walker_ESCC_TopMarker.csv')

write_rds(harmonySobj, 'data/walker_ESCC.rds')

# walker dataset read in
load('data/walker_OAC.DGE.RData')

OAC.DGE[1:5, 1:5]

CreateSeuratObject(OAC.DGE,
                   names.delim = '_',
                   min.cells = 10,
                   min.features = 3) -> sobj

sobj@assays$RNA[1:5,1:5]

sample <- read_tsv('data/walker_sample.txt')

sample_idents <- sample$sample

sobj <- subset(sobj, idents = sample_idents)

VlnPlot(sobj, 'nCount_RNA', log = TRUE)

sobj <- subset(sobj, nCount_RNA > 200)

metadata <- sobj@meta.data

metadata$genotype <- 'II'

metadata$genotype[which(metadata$orig.ident %in% c('OAC1411T', 'OAC174T'))] <- 'IT'

metadata$genotype[which(metadata$orig.ident == 'OAC132T')] <- 'TT'

ggplot(metadata)+
  geom_bar(aes(x = genotype, fill = genotype))

write_csv(metadata, 'results/walker_meta.csv')

sobj@meta.data <- metadata

SetIdent(sobj, value = 'genotype') -> sobj

load('data/walker_Cell.Metrics.RData')
write.csv(Metrics, 'results/walker_pub_meta.csv')
