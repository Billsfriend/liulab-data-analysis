# 2022.3.24

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load cell cycle markers
load("data/cycle.rda")

# load seurat objects
tpm10x <- readRDS("data/cell/tpm10x.rds")
tpmSmrt <- readRDS("data/cell/tpmSmart.rds")

# merge across platform
tpmMerge <- merge(tpm10x, tpmSmrt)

tpmMerge <- CellCycleScoring(tpmMerge, s.features = s_genes, g2m.features = g2m_genes)

# NormalTpm <- NormalizeData(NatureTpm)
# tpmMerge data need not normalize by counts per million but only log
tpmMerge[["RNA"]] <- CreateAssayObject(counts = log(tpmMerge[["RNA"]]@counts+1))
head(tpmMerge@assays$RNA@counts)

# Identify the most variable 2000 genes
tpmMerge <- FindVariableFeatures(tpmMerge, 
                                selection.method = "vst",
                                nfeatures = 2000, 
                                verbose = FALSE)

# Scale the counts
tpmMerge <- ScaleData(tpmMerge)

# Perform PCA
tpmMerge <- RunPCA(tpmMerge)

# Plot the PCA colored by I232T genotype
tpmMerge <- SetIdent(tpmMerge, value = 'ITgeno')
DimPlot(tpmMerge,
        reduction = "pca")

# Plot the PCA colored by sequence platforms
tpmMerge <- SetIdent(tpmMerge, value = 'Platform')
DimPlot(integMerge,
        reduction = "pca")

tpmMerge@meta.data$factor <- NA
meta <- tpmMerge@meta.data
meta$factor[which(meta$Platform=='10X'&meta$ITgeno=='IT')] <- 'IT10x'
meta$factor[which(meta$Platform=='10X'&meta$ITgeno=='II')] <- 'II10x'
meta$factor[which(meta$Platform=='Smart-seq2'&meta$ITgeno=='IT')] <- 'ITsmart'
meta$factor[which(meta$Platform=='Smart-seq2'&meta$ITgeno=='II')] <- 'IIsmart'
tpmMerge@meta.data <- meta

# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_tenx <- SplitObject(tpmMerge, split.by = "factor")
write_rds(split_tenx,'data/cell/split4merge.rds')

for (i in 1:length(split_tenx)) {
  #split_tenx[[i]] <- CellCycleScoring(split_tenx[[i]], g2m.features=g2m_genes, s.features=s_genes)
  split_tenx[[i]] <- SCTransform(split_tenx[[i]])
}

# Check which assays are stored in objects
split_tenx$IT10x@assays

# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_tenx, nfeatures = 3000)

# Prepare the SCT list object for integration
split_tenx <- PrepSCTIntegration(object.list = split_tenx, anchor.features = integ_features)

# Find best buddies - can take a while to run, ~42m28s
# progression will stay in 0 % while running, don't panic!
integ_anchors <- FindIntegrationAnchors(object.list = split_tenx, normalization.method = "SCT", anchor.features = integ_features)

# Integrate across conditions. can cost big space
integMerge <- IntegrateData(anchorset = integ_anchors, normalization.method = "SCT")

# save integrated file
write_rds(integMerge, 'data/cell/integMerge.rds')
