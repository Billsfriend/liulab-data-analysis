# 2022.6.13

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(org.Hs.eg.db)

# Load cell cycle markers
load("data/cycle.rda")

# load seurat objects
Srt <- readRDS("data/HolmesGC.rds")

Srt <- CellCycleScoring(Srt, s.features = s_genes, g2m.features = g2m_genes)

Srt <- NormalizeData(Srt)

# Identify the most variable 2000 genes
Srt <- FindVariableFeatures(Srt)

# Scale the counts
Srt <- ScaleData(Srt)

# Perform PCA
Srt <- RunPCA(Srt)

# Plot the PCA colored by cell cycle phase
DimPlot(Srt,
        reduction = "pca",
        split.by = 'Phase')

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_smart <- SplitObject(Srt, split.by = "Phase")

for (i in 1:length(split_smart)) {
  split_smart[[i]] <- CellCycleScoring(split_smart[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_smart[[i]] <- SCTransform(split_smart[[i]])
}

# Check which assays are stored in objects
split_smart$G2M@assays

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_smart, nfeatures = 3000)

# Prepare the SCT list object for integration
split_smart <- PrepSCTIntegration(object.list = split_smart, 
                                   anchor.features = integ_features)

# Find best buddies - can take a while to run, 26s
# progression will stay in 0 % while running, don't panic!
integ_anchors <- FindIntegrationAnchors(object.list = split_smart,
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)

# Integrate across conditions. can cost big space
smart_integrated <- IntegrateData(anchorset = integ_anchors, 
                                   normalization.method = "SCT")

# save integrated file
write_rds(smart_integrated, 'data/HolmesGC.rds')
