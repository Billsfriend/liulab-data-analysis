# devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
install.packages('harmony')
library(harmony)
library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(data.table)

# our data from Alevin
data <- fread('data/matrixSeurat.tsv.gz')
data <- data[2:60280]
data[1:5, 1:5]

barcodes <- fread('data/read10x/barcodes.tsv.gz', header = FALSE)

colnames(data) <- barcodes$V1

features <- fread('data/read10x/features.tsv.gz', header = FALSE)

rownames(data) <- features$V1

SLE_PBMC_Alevin <- CreateSeuratObject(data)

head(SLE_PBMC_Alevin@assays$RNA@data)


# quality control of auto-generated metadata
# calculate ratio of mitochondria genes
SLE_PBMC_Alevin$mitoRatio <- PercentageFeatureSet(object = SLE_PBMC_Alevin, pattern = "^MT-")
SLE_PBMC_Alevin$mitoRatio <- SLE_PBMC_Alevin@meta.data$mitoRatio / 100

# extract auto-generated data
orimeta <- SLE_PBMC_Alevin@meta.data
orimeta$log10GenesPerUMI <- log10(orimeta$nFeature_RNA) / log10(orimeta$nCount_RNA)
orimeta <- orimeta %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# delete column by dplyr
orimeta <- dplyr::select(orimeta, -orig.ident)
orimeta$barcodes <- rownames(orimeta)
# merge auto-generate and attached metadata
SLE_PBMC_Alevin@meta.data <- orimeta


# Visualize and quality control ---------
VlnPlot(SLE_PBMC_Alevin,
        features = 'nUMI') +
  ggtitle("Number of UMI per cell")

VlnPlot(SLE_PBMC_Alevin,
        features = 'nGene') +
  ggtitle("Number of genes per cell")

VlnPlot(SLE_PBMC_Alevin,
        features = 'mitoRatio') +
  ggtitle("ratio of mitochondria gene per cell")

SLE_PBMC_Alevin <- subset(SLE_PBMC_Alevin, mitoRatio < 0.10 & nGene > 200 & nGene < 2500)

# save data
write_rds(SLE_PBMC_Alevin, 'data/SLE_PBMC_Alevin.rds')

# processed data from publication ----------
mtx <- ReadMtx(mtx = 'data/read10x/matrix.mtx',
               cells = 'data/read10x/barcodes.tsv',
               features = 'data/read10x/features.tsv')

pub_srt <- CreateSeuratObject(mtx)

VlnPlot(pub_srt, features = 'nCount_RNA')

VlnPlot(pub_srt, features = 'nFeature_RNA')

# calculate ratio of mitochondria genes
pub_srt$mitoRatio <- PercentageFeatureSet(object = pub_srt, pattern = "^MT-")
pub_srt$mitoRatio <- pub_srt@meta.data$mitoRatio / 100

VlnPlot(pub_srt, features = 'mitoRatio')

pub_srt <- subset(pub_srt, mitoRatio < 0.10 & nFeature_RNA > 200 & nFeature_RNA < 2500)

write_rds(pub_srt, 'data/SLE_PBMC_Pub.rds')

# PCA SCT -----------
load("../CRC-I/data/cycle.rda")

pub_srt %>%
  CellCycleScoring(s.features = s_genes,
                   g2m.features = g2m_genes,
                   set.ident = TRUE) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA() ->
  pub_srt

DimPlot(pub_srt)

pub_srt %>%
  SplitObject() %>%
  lapply(SCTransform) ->
  split_pub

split_pub %>%
  SelectIntegrationFeatures(nfeatures = 3000) ->
  integ_features

split_pub %>%
  PrepSCTIntegration(anchor.features = integ_features) %>%
  FindIntegrationAnchors(
    normalization.method = 'SCT',
    anchor.features = integ_features
  ) %>%
  IntegrateData(normalization.method = 'SCT') %>%
  RunPCA() ->
  integrated_pub

# Alevin output
SLE_PBMC_Alevin %>%
  CellCycleScoring(s.features = s_genes,
                   g2m.features = g2m_genes,
                   set.ident = TRUE) %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA() ->
  SLE_PBMC_Alevin

DimPlot(SLE_PBMC_Alevin)

SLE_PBMC_Alevin %>%
  SplitObject() %>%
  lapply(SCTransform) ->
  split_alevin

split_alevin %>%
  SelectIntegrationFeatures(nfeatures = 3000) ->
  integ_features2

split_alevin %>%
  PrepSCTIntegration(anchor.features = integ_features2) %>%
  FindIntegrationAnchors(
    normalization.method = 'SCT',
    anchor.features = integ_features2
  ) %>%
  IntegrateData(normalization.method = 'SCT') %>%
  RunPCA() ->
  integrated_alevin

# identify markers --------
DimPlot(integrated_pub)

DimPlot(integrated_alevin)

integrated_pub %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunTSNE() ->
  integrated_pub

DimPlot(integrated_pub)

integrated_alevin %>%
  FindNeighbors() %>%
  FindClusters() %>%
  RunTSNE() ->
  integrated_alevin

DimPlot(integrated_alevin)
