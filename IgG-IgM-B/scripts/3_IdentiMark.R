# Single-cell RNA-seq analysis - clustering analysis

# BiocManager::install('multtest')
# `bash sudo apt install libfftw3-dev`
# install.packages('metap')
# Load libraries
library(Seurat)
library(tidyverse)
library(cowplot)
library(data.table)
library(scales)

# stop using view() in GUI for big matrix, it is VERY SLOW
# Load integrated data
Srt <- readRDS("data/HolmesLZ.rds")

# Run PCA
Srt <- RunPCA(object = Srt)

# Plot PCA
PCAPlot(Srt,
        split.by = "ident")  

# Run UMAP with 40 PCs
Srt <- RunUMAP(Srt,
               dims = 1:40,
               reduction = "pca")

# Plot UMAP                             
DimPlot(Srt)

DimPlot(Srt,
        split.by = "Phase")

# Explore heatmap of PCs
DimHeatmap(Srt, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = Srt[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = Srt, 
          ndims = 40)

# Determine the K-nearest neighbor graph
Srt <- FindNeighbors(object = Srt,                                   dims = 1:40)

# Determine the clusters for various resolutions                  
Srt <- FindClusters(object = Srt,
                    resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Assign identity of cluster by res 0.8
Idents(object = Srt) <- "integrated_snn_res.0.8"

# Plot the UMAP at res 0.8
DimPlot(Srt,
        reduction = "umap",
        label = TRUE,
        label.size = 6) +
  scale_color_viridis_d(option = 'turbo', begin = 0.1)+
  labs(title = 'UMAP 0.8')

# Assign identity of clusters at resolution of 0.6
Idents(object = Srt) <- "integrated_snn_res.0.6"

# Plot the UMAP 0.6
DimPlot(Srt,
        reduction = "umap",
        label = TRUE,
        label.size = 6)+ scale_color_viridis_d(option = 'turbo', begin = 0.1)+
  labs(title = 'UMAP 0.6')

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(Srt, 
                     vars = c("ident", "Phase")) %>%
  dplyr::count(ident, Phase) %>%
  tidyr::spread(ident, n)

# View table
# View(n_cells)

# Explore whether clusters segregate by cell cycle phase
DimPlot(Srt,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in Srt@meta.data
metrics <-  c("nUMI", "nGene", "S.Score", "G2M.Score")

FeaturePlot(Srt, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)


# Examine PCA results 
print(Srt[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(object = Srt, 
        reduction = "umap", 
        label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(Srt) <- "RNA"

# Normalize RNA data for visualization purposes
Srt <- NormalizeData(Srt)

# Plasmacytoid dendritic cell markers
# FeaturePlot(Srt, 
#             reduction = "umap", 
#             features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"), 
#             order = TRUE,
#             min.cutoff = 'q10', 
#             label = TRUE)

# explicitly set our default assay
# we want to use the original counts and not the integrated data
DefaultAssay(Srt) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(Srt,
                                                   ident.1 = 0,
                                                   grouping.var = "Phase",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
# Since each cell is being treated as a replicate
# this will result in inflated p-values within each group!

annotations <- read.csv("data/annotation.csv")

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(Srt,
                       ident.1 = cluster,
                       grouping.var = "Phase",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:11), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (G1_avg_log2FC + S_avg_log2FC + G2M_avg_log2FC) /3) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)

# change ident to see count-base classification
Srt <- SetIdent(Srt, value = 'exClass')

Srt <- ScaleData(Srt)

# Plot interesting marker gene expression
FeaturePlot(object = Srt, 
            features = c("IGHG1",'IGHM'))

DotPlot(object = Srt, 
          features = c("IGHG1",'IGHM'))

VlnPlot(object = Srt, 
        features = c("IGHG1",'IGHM'))

DotPlot(object = Srt, 
          features = c("LY9",'CD86','ICOSLG','CD40','CD9','ICAM2'))

# markers in Holmes2020 paper
DotPlot(object = Srt, 
        features = c("CXCR4",'CD83','PRDM1','CCR6'))

FeaturePlot(object = Srt, 
        features = c("CXCR4",'CD83','PRDM1','CCR6','CD9'))

VlnPlot(Srt, features = 'CD40')

# run tSNE on PCA reduction
Srt <- RunTSNE(Srt)

# extract high express IgG1/IgM clusters
refine <- subset(Srt, idents = c('1','2','LZ2','7','LZ1'))

library(tidyseurat)

refine$class <- NA

# rename for class
refine$class[which(refine$integrated_snn_res.0.8 %in% c('1','2','5'))] <- 'IgM'

refine$class[which(refine$integrated_snn_res.0.8 %in% c('0','7'))] <- 'IgG1'

metadata <- refine@meta.data

# examine in interactive table
library(DT)
datatable(metadata)

# set as ident
refine <- SetIdent(refine, value = 'class')

# save
write_rds(refine, 'data/LZclass.rds')

# save labeling seurat
write_rds(Srt, 'data/HolmesLZ.rds')

# save top10 markers
write_csv(top10, 'results/Cluster-HolmesGC-10clustersMarker.csv')
