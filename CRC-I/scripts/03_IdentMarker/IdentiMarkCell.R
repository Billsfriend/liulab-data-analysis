# Single-cell RNA-seq analysis - clustering analysis

# Load libraries
library(Seurat)
library(data.table)
library(tidyverse)

# stop using view() in GUI for big matrix, it is VERY SLOW
# Load integrated data
sobj <- readRDS("data/cell/tpm_smart.rds")

# Run UMAP
# default need install python package umap-learn
sobj <- RunUMAP(sobj,
  dims = 1:40,
)

# Plot UMAP
DimPlot(sobj)

DimPlot(sobj,
  split.by = "Phase"
)

# Plot the elbow plot
ElbowPlot(
  object = sobj,
  ndims = 40
)

# Determine the K-nearest neighbor graph
sobj <- FindNeighbors(
  object = sobj,
  dims = 1:40
)

# Determine the clusters for various resolutions
sobj <- FindClusters(
  object = sobj,
  resolution = 0.6
)

# Plot the UMAP
DimPlot(sobj,
  label = TRUE,
  label.size = 6
) + scale_color_viridis_d(option = "turbo", begin = 0.1)


# UMAP of cells in each cluster by genotypes
DimPlot(sobj,
  label = TRUE,
  split.by = "ITgeno"
) + NoLegend()

# Explore whether clusters segregate by cell cycle phase
DimPlot(sobj,
  label = TRUE,
  split.by = "Phase"
) + NoLegend()

# Determine metrics to plot present in sobj@meta.data
metrics <- c("nUMI", "nGene", "S.Score", "G2M.Score")

FeaturePlot(sobj,
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
pc_data <- FetchData(sobj,
  vars = columns
)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(sobj,
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
print(sobj[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(
  object = sobj,
  reduction = "umap",
  label = TRUE
) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(sobj) <- "RNA"

# Normalize RNA data for visualization purposes
sobj <- NormalizeData(sobj, verbose = FALSE)

# Plasmacytoid dendritic cell markers
# FeaturePlot(sobj,
#             reduction = "umap",
#             features = c("IL3RA", "GZMB", "SERPINF1", "ITM2C"),
#             order = TRUE,
#             min.cutoff = 'q10',
#             label = TRUE)

# explicitly set our default assay
# we want to use the original counts and not the integrated data
DefaultAssay(sobj) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(sobj,
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
  FindConservedMarkers(sobj,
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
conserved_markers <- map_dfr(c(0:21), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>%
  mutate(avg_fc = (II_avg_log2FC + IT_avg_log2FC) / 2) %>%
  group_by(cluster_id) %>%
  top_n(
    n = 10,
    wt = avg_fc
  )

# Visualize top 10 markers per cluster
View(top10)

# Plot interesting marker gene expression
FeaturePlot(
  object = sobj,
  features = c("DNAJA1", "DNAJB1", "HSPA8"),
  order = TRUE,
  min.cutoff = "q10",
  split.by = "ITgeno"
)

# Vln plot - cluster
VlnPlot(
  object = sobj,
  features = c("HLA-G", "HLA-J", "HLA-H")
)

RidgePlot(
  object = sobj,
  features = c("HLA-G", "HLA-J", "HLA-H"), group.by = "Tissue"
)

VlnPlot(
  object = sobj,
  features = c("DNAJA1", "DNAJB1", "HSPA8")
)

RidgePlot(
  object = sobj,
  features = c("DNAJA1", "DNAJB1", "HSPA8")
)

# save labeling seurat
write_rds(sobj, "data/cell/IdentMerge.rds")

splitMerge <- SplitObject(sobj, split.by = "Tissue")
write_rds(splitMerge[["T"]], "data/cell/mergeTumor.rds")
write_rds(splitMerge[["P"]], "data/cell/mergeBlood.rds")
write_rds(splitMerge[["N"]], "data/cell/mergeNormal.rds")

# save top10 markers
write_csv(top10, "results/cluster/CellMerge-22clustersMarker.csv")

# show difference between our clusters and Zhang's clusters
sobj@meta.data %>%
  ggplot(aes(x = Global_Cluster, fill = integrated_snn_res.0.4)) +
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
  ggtitle("NCells") +
  scale_fill_viridis_d()


DimPlot(sobj, split.by = "ITgeno", group.by = "majorCluster", pt.size = 1) + scale_color_viridis_d(option = "turbo", begin = 0.1)
