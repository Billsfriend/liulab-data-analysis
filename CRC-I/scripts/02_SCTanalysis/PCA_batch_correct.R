# 2022.3.24
# edit in 2022.9.29

# Load libraries
library(Seurat)
library(harmony)
library(tidyverse)

# load seurat objects
sobj <- readRDS("data/cell/tpmSmart.rds")

# sobj data need not normalize by counts per million but only log
sobj[["RNA"]] <- CreateAssayObject(counts = log(sobj[["RNA"]]@counts + 1))
sobj@assays$RNA@counts[1:5, 1:5]
sobj@assays$RNA@data[1:5, 1:5]

# Load cell cycle markers
Seurat::cc.genes.updated.2019 -> cc.genes

sobj %>%
  CellCycleScoring(
    cc.genes$g2m.genes,
    cc.genes$s.genes
  ) -> sobj

# Identify the most variable 2000 genes
sobj <- FindVariableFeatures(sobj, verbose = FALSE)

# Scale the counts
sobj <- ScaleData(sobj)
sobj <- SetIdent(sobj, value = "genotype")

# Perform PCA
sobj <- RunPCA(sobj)

# Plot the PCA colored by cell cycle phase
DimPlot(sobj,
  reduction = "pca"
)

RunHarmony(sobj, 'Sample') -> sobj

DimPlot(sobj, reduction = 'harmony')

# save integrated file
write_rds(sobj, "data/cell/tpm_smart.rds")
