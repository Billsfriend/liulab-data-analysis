library(data.table)
library(Seurat)
library(tidyverse)
library(harmony)

# create Seurat as expr-matrix -------
data <- fread('data/Wang2022/GSE110009_metastatic_Colon_TPM.txt.gz')

data[,-1][1:5,1:5]

data <- unique(data, by = "V1")

rownames(data) <- data$V1

sobj <- CreateSeuratObject(data[, -1],
                           row.names = rownames(data))

head(sobj)

sobj@assays$RNA@counts[1:5, 1:5]

# auto-assign ident based on barcode naming
sobj@active.ident

sobj %>%
  PercentageFeatureSet(pattern = "^MT-") ->
  sobj$mitoRatio

sobj %>%
  PercentageFeatureSet(pattern = "^HLA-") ->
  sobj$hlaRatio

VlnPlot(sobj, 'nFeature_RNA')
VlnPlot(sobj, 'hlaRatio')

sobj <- subset(sobj, nFeature_RNA > 500)

sobj@meta.data -> metadata
metadata$barcode <- rownames(metadata)
metadata <- separate(metadata,
              col = barcode,
              sep = '_',
              into = c('patient', 'tissue', 'type', 'id'))

metadata <- unite(metadata,
                  patient,
                  tissue,
                  type,
                  col = 'batch',
                  sep = '_',
                  remove = FALSE)

metadata$genotype <- NA

metadata$genotype[which(metadata$patient == 'P5')] <- 'TT' 

metadata$genotype[which(metadata$patient %in% c('P3','P4','P7','P9'))] <- 'IT' 

metadata$genotype[which(metadata$patient %in% c('P1','P2','P6','P8'))] <- 'II' 

# visualize cell number
metadata %>%
  ggplot() +
  geom_bar(aes(x = genotype, fill = genotype)) +
  labs(title = 'Tumor cell number')

sobj@meta.data <- metadata

write_rds(sobj, 'data/Wang2022/wang2022.rds')

sobj <- read_rds('data/Wang2022/wang2022.rds')

# retain only tumor cells
tumorBarcode <- str_which(sobj$tissue, 'T')

sobj <- subset(sobj, cells = tumorBarcode)

write_rds(sobj, 'data/Wang2022/wangTumor.rds')

# normalize TPM to logTPM+1
sobj[["RNA"]] <- CreateAssayObject(counts = log(sobj[["RNA"]]@counts+1))

# essential commands
sobj %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() ->
  sobj

DimPlot(sobj, group.by = 'genotype')

# try use harmony to clear batch effect ------
harmonySobj <- RunHarmony(sobj, 
                          'genotype')

DimPlot(harmonySobj,
        reduction = 'harmony',
        group.by = 'genotype')

# run UMAP and clustering
harmonySobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity() ->
  harmonySobj

DimPlot(harmonySobj, reduction = 'umap')
DimPlot(harmonySobj, 
        reduction = 'umap',
        split.by = 'genotype')

RunTSNE(harmonySobj) -> harmonySobj
TSNEPlot(harmonySobj)

# get conserved markers for any given cluster ----------
annotations <- read.csv("data/annotation.csv")

# Iterate function across desired clusters
allMarkers <- FindAllMarkers(harmonySobj,
                             only.pos = TRUE,
                             min.pct = 0.25)

allMarkers %>%
  group_by(cluster) %>%
  top_n(n = -10, wt = p_val_adj) ->
  topMarkers

write_csv(topMarkers, 'results/cluster/wang2022allTopMarker.csv')