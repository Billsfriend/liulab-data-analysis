# Single-cell RNA-seq analysis - clustering analysis

# BiocManager::install('multtest')
# `bash sudo apt install libfftw3-dev`
# install.packages('metap')
# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)
library(data.table)
library(scales)

# stop using view() in GUI for big matrix, it is VERY SLOW
# Load integrated data
integNat <- readRDS("data/nature/integNat.rds")

# Run PCA
integNat <- RunPCA(object = integNat)

# Plot PCA
PCAPlot(integNat,
  split.by = "ITgeno"
)

# Run UMAP
integNat <- RunUMAP(integNat,
  dims = 1:40,
  reduction = "pca"
)

# Plot UMAP
DimPlot(integNat)

DimPlot(integNat,
  split.by = "ITgeno"
)

# Explore heatmap of PCs
DimHeatmap(integNat,
  dims = 1:9,
  cells = 500,
  balanced = TRUE
)

# Printing out the most variable genes driving PCs
print(
  x = integNat[["pca"]],
  dims = 1:10,
  nfeatures = 5
)

# Plot the elbow plot
ElbowPlot(
  object = integNat,
  ndims = 40
)

# Determine the K-nearest neighbor graph
integNat <- FindNeighbors(
  object = integNat,
  dims = 1:40
)

# Determine the clusters for various resolutions
integNat <- FindClusters(
  object = integNat,
  resolution = c(0.4, 0.6, 0.8, 1.0, 1.4)
)

# Assign identity of clusters
Idents(object = integNat) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(integNat,
  reduction = "umap",
  label = TRUE,
  label.size = 6
) + scale_color_viridis_d(option = "turbo", begin = 0.1)

# Assign identity of clusters at resolution of 0.4
Idents(object = integNat) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(integNat,
  reduction = "umap",
  label = TRUE,
  label.size = 6
) + scale_color_viridis_d(option = "turbo", begin = 0.1)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(integNat,
  vars = c("ident", "ITgeno")
) %>%
  dplyr::count(ident, ITgeno) %>%
  tidyr::spread(ident, n)

# View table
# View(n_cells)

# UMAP of cells in each cluster by genotypes
DimPlot(integNat,
  label = TRUE,
  split.by = "ITgeno"
) + NoLegend()

# Explore whether clusters segregate by cell cycle phase
DimPlot(integNat,
  label = TRUE,
  split.by = "Phase"
) + NoLegend()

# Determine metrics to plot present in integNat@meta.data
metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score")

FeaturePlot(integNat,
  reduction = "umap",
  features = metrics,
  pt.size = 0.4,
  order = TRUE,
  min.cutoff = "q10",
  label = TRUE
)

# Defining the information in the seurat object of interest
columns <- c(
  paste0("PC_", 1:16),
  "ident",
  "UMAP_1", "UMAP_2"
)

# Extracting this data from the seurat object
pc_data <- FetchData(integNat,
  vars = columns
)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(integNat,
  vars = c("ident", "UMAP_1", "UMAP_2")
) %>%
  group_by(ident) %>%
  summarise(x = mean(UMAP_1), y = mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
# map(paste0("PC_", 1:16), function(pc){
#   ggplot(pc_data,
#          aes(UMAP_1, UMAP_2)) +
#     geom_point(aes_string(color=pc),
#                alpha = 0.7) +
#     scale_color_gradient(guide = "none",
#                          low = "grey90",
#                          high = "blue")  +
#     geom_text(data=umap_label,
#               aes(label=ident, x, y)) +
#     ggtitle(pc)
# }) %>%
#   plot_grid(plotlist = .)

# Examine PCA results
print(integNat[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(
  object = integNat,
  reduction = "umap",
  label = TRUE
) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(integNat) <- "RNA"

# Normalize RNA data for visualization purposes
integNat <- NormalizeData(integNat, verbose = FALSE)

# Plasmacytoid dendritic cell markers
# FeaturePlot(integNat,
#             reduction = "umap",
#             features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"),
#             order = TRUE,
#             min.cutoff = 'q10',
#             label = TRUE)

# explicitly set our default assay
# we want to use the original counts and not the integrated data
DefaultAssay(integNat) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(integNat,
  ident.1 = 0,
  grouping.var = "ITgeno",
  only.pos = TRUE,
  logfc.threshold = 0.25
)
# Since each cell is being treated as a replicate
# this will result in inflated p-values within each group!

annotations <- read.csv("data/annotation.csv")

# Combine markers with gene descriptions
cluster0_ann_markers <- cluster0_conserved_markers %>%
  rownames_to_column(var = "gene") %>%
  left_join(
    y = unique(annotations[, c("gene_name", "description")]),
    by = c("gene" = "gene_name")
  )

View(cluster0_ann_markers)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster) {
  FindConservedMarkers(integNat,
    ident.1 = cluster,
    grouping.var = "ITgeno",
    only.pos = TRUE
  ) %>%
    rownames_to_column(var = "gene") %>%
    left_join(
      y = unique(annotations[, c("gene_name", "description")]),
      by = c("gene" = "gene_name")
    ) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:12), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>%
  mutate(avg_fc = (II_avg_log2FC + IT_avg_log2FC + TT_avg_log2FC) / 3) %>%
  group_by(cluster_id) %>%
  top_n(
    n = 10,
    wt = avg_fc
  )

# Visualize top 10 markers per cluster
View(top10)

# Plot interesting marker gene expression
FeaturePlot(
  object = labelNat,
  features = c("DNAJA1", "DNAJB1", "HSPA8"),
  order = TRUE,
  min.cutoff = "q10",
  split.by = "ITgeno"
)

# Vln plot - cluster
VlnPlot(
  object = labelNat,
  features = c("HLA-G", "HLA-J", "HLA-H")
)

RidgePlot(
  object = labelNat,
  features = c("HLA-G", "HLA-J", "HLA-H"), group.by = "tissue"
)

VlnPlot(
  object = labelNat,
  features = c("DNAJA1", "DNAJB1", "HSPA8")
)

RidgePlot(
  object = labelNat,
  features = c("DNAJA1", "DNAJB1", "HSPA8")
)

# run tSNE on PCA reduction
integNat <- RunTSNE(integNat)

DimPlot(CountNat, reduction = "tsne", group.by = "SingleR.MI.labels")

# save labeling seurat
write_rds(integNat, "data/nature/CountNat.rds")

# save top10 markers
write_csv(top10, "results/cluster/NatCount-13clustersMarker.csv")

# show difference between our clusters and Zhang's clusters
integNat@meta.data %>%
  ggplot(aes(x = majorCluster, fill = integrated_snn_res.0.8)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells") +
  scale_fill_viridis_d()


DimPlot(labelNat, split.by = "ITgeno", group.by = "majorCluster", pt.size = 1) + scale_color_viridis_d(option = "turbo", begin = 0.1)
