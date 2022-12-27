library(Seurat)
library(tidyverse)
library(slingshot)
library(DropletUtils)

labelNat <- read_rds('data/nature/CountNat.rds')

DimPlot(labelNat, group.by = 'majorCluster')+
  scale_color_viridis_d(option = 'turbo', begin = 0.1)

DimPlot(labelNat)+
  scale_color_viridis_d(option = 'turbo', begin = 0.1)


byI232T <- SplitObject(labelNat, split.by = 'ITgeno')


# separately write 10x counts directories for each genotype
write10xCounts(path = 'data/nature_CRC_Tcell_II
',
x = byI232T[["II"]]@assays[["RNA"]]@counts)

write10xCounts(path = 'data/nature_CRC_Tcell_IT
',
x = byI232T[["IT"]]@assays[["RNA"]]@counts)

write10xCounts(path = 'data/nature_CRC_Tcell_TT
',
x = byI232T[["TT"]]@assays[["RNA"]]@counts)

# majorCluster <- labelNat@meta.data$majorCluster
# umap <- labelNat@reductions$umap
# 
# scExperNat <- as.SingleCellExperiment(labelNat)
# 
# reducedDim(scExperNat)
# 
# scExperNat <- slingshot(scExperNat, clusterLabels = 'majorCluster', reducedDim= 'PCA')
