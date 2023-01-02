library(Seurat)
library(ggpubr)
library(tidyverse)

sobj <- read_rds('IgG-IgM-B/data/Liu-Qi_2022_LLPC/aggr.rds')

sobj@meta.data <- sobj@meta.data %>%
  rownames_to_column('rownames') %>%
  separate(orig.ident, into = c('time', 'tissue'), remove = FALSE) %>%
  mutate(tissue = case_when(
    str_detect(orig.ident, 'bm') ~ 'bone_marrow',
    TRUE ~ 'spleen'
  )) %>%
  column_to_rownames('rownames')

t_features <- c('Cd3d','Cd3e','Cd3g','Cd247','Trac','Trbc1','Trbc2','Trdc')

sobj <- sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA()

DimPlot(sobj)

sobj <- harmony::RunHarmony(sobj, 'orig.ident')

DimPlot(sobj, reduction = 'harmony')

sobj <- sobj %>%
  RunUMAP(dims = 1:20, reduction = 'harmony') %>%
  FindNeighbors(dims = 1:20, reduction = 'harmony') %>%
  FindClusters()

DimPlot(sobj)

rownames(sobj) %>%
  str_subset('Ighg')

VlnPlot(sobj, c('Ighg1', 'Ighm', 'Igha'), pt.size = 0)

DotPlot(sobj, features = c('Ighg1', 'Ighg2b', 'Ighg2c', 'Ighg3', 'Ighm', 'Igha'), scale = FALSE)

Idents(sobj) <- 'tissue'

DotPlot(sobj,
        features = c('Prdm1', 'Sdc1'),
        scale.by = 'size',
        scale = FALSE,
        group.by = 'time',
        idents = 'bone_marrow') +
  labs(title = 'Blimp & Cd138 expression in bone marrow plasma cells',
       x = 'Genes',
       y = 'Time')

bm <- last_plot()

DotPlot(sobj,
        features = c('Prdm1', 'Sdc1'),
        scale.by = 'size',
        scale.min = 10,
        scale = FALSE,
        group.by = 'time',
        idents = 'spleen') +
  labs(title = 'Blimp & Cd138 expression in spleen plasma cells',
       x = 'Genes',
       y = 'Time')

sp <- last_plot()

bm + NoLegend() + sp

DotPlot(sobj,
        features = t_features,
        scale.by = 'size',
        scale.max = 2,
        scale = FALSE,
        group.by = 'time',
        idents = 'bone_marrow') +
  labs(title = 'TCR expression in bone marrow plasma cells',
       x = 'Genes',
       y = 'Time')

bm <- last_plot()

DotPlot(sobj,
        features = t_features,
        scale.by = 'size',
        scale.max = 2,
        scale = FALSE,
        group.by = 'time',
        idents = 'spleen') +
  labs(title = 'TCR expression in spleen plasma cells',
       x = 'Genes',
       y = 'Time')

sp <- last_plot()

bm + NoLegend() + sp

sobj <- AddModuleScore(sobj, 
                       features = list(t_features),
                       name = 'TCR_score')
