# Single-cell RNA-seq analysis - clustering analysis

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load integrated data
integSmart <- readRDS("data/cell/integSmart.rds")
# integSmart <- SetIdent(integSmart, value = 'ITgeno')

# Run PCA
integSmart <- RunPCA(object = integSmart)

# Plot PCA
PCAPlot(integSmart,
        split.by = "ITgeno")  

# Run UMAP
integSmart <- RunUMAP(integSmart, 
                    dims = 1:40,
                    reduction = "pca")

# Plot UMAP                             
DimPlot(integSmart)

DimPlot(integSmart,
        split.by = "ITgeno")

# Explore heatmap of PCs
DimHeatmap(integSmart, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = integSmart[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = integSmart, 
          ndims = 40)

# Determine the K-nearest neighbor graph
integSmart <- FindNeighbors(object = integSmart, 
                                   dims = 1:40)

# Determine the clusters for various resolutions
integSmart <- FindClusters(object = integSmart, resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolutions
integSmart@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = integSmart) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(integSmart,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Assign identity of clusters at resolution of 0.4
Idents(object = integSmart) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(integSmart,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(integSmart, 
                     vars = c("ident", "ITgeno")) %>%
  dplyr::count(ident, ITgeno) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by genotypes
DimPlot(integSmart, 
        label = TRUE, 
        split.by = "ITgeno")  + NoLegend()

# Explore whether clusters segregate by cell cycle phase
DimPlot(integSmart,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in integSmart@meta.data
metrics <-  c("filter.nUMI", "filter.nGene", "S.Score", "G2M.Score")

FeaturePlot(integSmart, 
            reduction = "umap", 
            features = metrics,
            pt.size = 0.4, 
            order = TRUE,
            min.cutoff = 'q10',
            label = TRUE)

# Defining the information in the seurat object of interest
columns <- c(paste0("PC_", 1:16),
             "ident",
             "UMAP_1", "UMAP_2")

# Extracting this data from the seurat object
pc_data <- FetchData(integSmart, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(integSmart, 
                        vars = c("ident", "UMAP_1", "UMAP_2"))  %>%
  group_by(ident) %>%
  summarise(x=mean(UMAP_1), y=mean(UMAP_2))

# Plotting a UMAP plot for each of the PCs
map(paste0("PC_", 1:16), function(pc){
  ggplot(pc_data, 
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=pc), 
               alpha = 0.7) +
    scale_color_gradient(guide = "none", 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(pc)
}) %>% 
  plot_grid(plotlist = .)

# Examine PCA results 
print(integSmart[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(object = integSmart, 
        reduction = "umap", 
        label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(integSmart) <- "RNA"

# Normalize RNA data for visualization purposes
# integSmart <- NormalizeData(integSmart, verbose = FALSE)

# mark cell expressing specific genes e.g. CD14+ monocytes
FeaturePlot(integSmart, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# explicitly set our default assay
# we want to use the original counts and not the integrated data
DefaultAssay(integSmart) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(integSmart,
                                                   ident.1 = 0,
                                                   grouping.var = "ITgeno",
                                                   only.pos = TRUE,
                                                   logfc.threshold = 0.25)
# Since each cell is being treated as a replicate
# this will result in inflated p-values within each group!

annotations <- read.csv("data/annotation.csv")

# Combine markers with gene descriptions 
cluster0_ann_markers <- cluster0_conserved_markers %>% 
  rownames_to_column(var="gene") %>% 
  left_join(y = unique(annotations[, c("gene_name", "description")]),
            by = c("gene" = "gene_name"))

View(cluster0_ann_markers)

# Create function to get conserved markers for any given cluster
get_conserved <- function(cluster){
  FindConservedMarkers(integSmart,
                       ident.1 = cluster,
                       grouping.var = "ITgeno",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:15), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (II_avg_log2FC + IT_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Plot interesting marker gene expression
FeaturePlot(object = integSmart, 
            features = c("FCGR2B"),
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

# Vln plot - cluster 
VlnPlot(object = integSmart, 
        features = c("FCGR2B"))

integSmart <- SetIdent(integSmart, value = 'ITgeno')

# save labeling seurat
write_rds(integSmart, 'data/cell/labelSmart.rds')

# save top10 markers
write_csv(top10, 'results/cluster/Smart-clustersMarker.csv')

# show difference between our clusters and Zhang's clusters
integSmart@meta.data %>% 
  ggplot(aes(x=Sub_Cluster, fill=integrated_snn_res.1)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

