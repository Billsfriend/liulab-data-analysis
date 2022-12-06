library(Seurat)
library(tidyverse)
library(tidyseurat)
library(SingleR)
library(ggpubr)

list.files('FPP_TRPV3/data/Dunlap-Rao_2022_CLE/', full.names = TRUE) -> f

map(f, Read10X) -> data

map(data, CreateSeuratObject) -> sobjlist

add_label <- function(x, label){
  x$sample <- label
  x}

str_extract(f, '(?<=//).+(?=_co)') -> samples

map2(sobjlist, samples, add_label) -> labellist

reduce(labellist, merge) -> sobj

write_rds(sobj, 'FPP_TRPV3/data/Dunlap-Rao_2022_CLE/sobj.rds', compress = 'gz')

# save to h5 and convert to anndata to analyze on FASTgenomics
SeuratDisk::SaveH5Seurat(sobj, 'FPP_TRPV3/data/dunlap_cle.h5seurat')

SeuratDisk::Convert('FPP_TRPV3/data/dunlap_cle.h5seurat', dest = 'h5ad')

# after merging ----------
read_rds('FPP_TRPV3/data/Dunlap-Rao_2022_CLE/sobj.rds') -> sobj

sobj %>%
  mutate(orig.ident = case_when(
    str_detect(sample, 'Healthy') ~ 'Healthy',
    str_detect(sample, 'Non') ~ 'Lupus_nonlesion',
    TRUE ~ 'Lupus_lesion'
  )) -> sobj

Idents(sobj) <- sobj$orig.ident

# pre-process -------------
sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() -> sobj

DimPlot(sobj)

harmony::RunHarmony(sobj, group.by.vars = 'sample') -> sobj

sobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity() ->
  sobj

DimPlot(sobj)

sobj <- subset(sobj, orig.ident != 'Lupus_nonlesion')

sample_n(sobj, 10000) -> small_sobj

per_clust_idented <- sobj %>%
  as.SingleCellExperiment() %>%
  SingleR::SingleR(
    ref = human_cell_ref,
    labels = human_cell_ref$label.main,
    clusters = sobj$seurat_clusters,
    de.method = "wilcox",
    BPPARAM = BiocParallel::MulticoreParam()
  )

table(per_clust_idented$pruned.labels)

plotScoreHeatmap(per_clust_idented,
                 show_colnames = TRUE,
                 show.pruned = TRUE)

match(sobj@meta.data$seurat_clusters,
      rownames(per_clust_idented)) -> match_index

sobj$singler_labels <- per_clust_idented$pruned.labels[match_index]

DimPlot(sobj, group.by = 'singler_labels')+
  ggtitle('Major cell types in CLE/HC skin')+
  scale_color_viridis_d(begin = 0.1, option = 'turbo')

DimPlot(sobj, group.by = 'fine_label')+
  ggtitle('Fine cell types in CLE/HC skin')+
  scale_color_viridis_d(begin = 0.1, option = 'turbo')

filter(sobj, singler_labels == 'Keratinocytes') %>%
  DimPlot(group.by = 'fine_label')+
  ggtitle('Fine cell types in CLE/HC keratinocytes')+
  scale_color_viridis_d(begin = 0.1, option = 'turbo')

ggsave('FPP_TRPV3/figures/umap.pdf')

sobj %>%
  unite('fine_label',
        c('singler_labels','seurat_clusters'),
        remove = FALSE) ->
  sobj

SetIdent(sobj, value = 'fine_label') -> sobj

write_rds(sobj, 'FPP_TRPV3/data/Dunlap-Rao_2022_CLE/sobj.rds', compress = 'gz')
