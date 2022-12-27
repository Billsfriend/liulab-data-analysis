# Single-cell RNA-seq analysis - clustering analysis

# Load libraries
library(Seurat)
library(tidyverse)
library(RCurl)
library(cowplot)

# Load integrated data
integ10x <- readRDS("data/cell/integ10x.rds")

# Run PCA
integ10x <- RunPCA(object = integ10x)

# Plot PCA
PCAPlot(integ10x,
        split.by = "ITgeno")  

# Run UMAP
integ10x <- RunUMAP(integ10x, 
                    dims = 1:40,
                    reduction = "pca")

# Plot UMAP                             
DimPlot(integ10x)

DimPlot(integ10x,
        split.by = "ITgeno")

# Explore heatmap of PCs
DimHeatmap(integ10x, 
           dims = 1:9, 
           cells = 500, 
           balanced = TRUE)

# Printing out the most variable genes driving PCs
print(x = integ10x[["pca"]], 
      dims = 1:10, 
      nfeatures = 5)

# Plot the elbow plot
ElbowPlot(object = integ10x, 
          ndims = 40)

# Determine the K-nearest neighbor graph
integ10x <- FindNeighbors(object = integ10x, 
                                   dims = 1:40)

# Determine the clusters for various resolutions            
integ10x <- FindClusters(object = integ10x,
                        resolution = c(0.4, 0.6, 0.8, 1.0, 1.4))

# Explore resolutions
integ10x@meta.data %>% 
  View()

# Assign identity of clusters
Idents(object = integ10x) <- "integrated_snn_res.0.8"

# Plot the UMAP
DimPlot(integ10x,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Assign identity of clusters at resolution of 0.4
Idents(object = integ10x) <- "integrated_snn_res.0.4"

# Plot the UMAP
DimPlot(integ10x,
        reduction = "umap",
        label = TRUE,
        label.size = 6)

# Extract identity and sample information from seurat object to determine the number of cells per cluster per sample
n_cells <- FetchData(integ10x, 
                     vars = c("ident", "ITgeno")) %>%
  dplyr::count(ident, ITgeno) %>%
  tidyr::spread(ident, n)

# View table
View(n_cells)

# UMAP of cells in each cluster by genotypes
DimPlot(integ10x, 
        label = TRUE, 
        split.by = "ITgeno")  + NoLegend()

# Explore whether clusters segregate by cell cycle phase
DimPlot(integ10x,
        label = TRUE, 
        split.by = "Phase")  + NoLegend()

# Determine metrics to plot present in integ10x@meta.data
metrics <-  c("filter.nUMI", "filter.nGene", "S.Score", "G2M.Score")

FeaturePlot(integ10x, 
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
pc_data <- FetchData(integ10x, 
                     vars = columns)

# Adding cluster label to center of cluster on UMAP
umap_label <- FetchData(integ10x, 
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
print(integ10x[["pca"]], dims = 1:5, nfeatures = 5)

DimPlot(object = integ10x, 
        reduction = "umap", 
        label = TRUE) + NoLegend()

# Select the RNA counts slot to be the default assay
DefaultAssay(integ10x) <- "RNA"

# Normalize RNA data for visualization purposes
# integ10x <- NormalizeData(integ10x, verbose = FALSE)

# mark CD14+ monocytes
FeaturePlot(integ10x, 
            reduction = "umap", 
            features = c("CD14", "LYZ"), 
            order = TRUE,
            min.cutoff = 'q10', 
            label = TRUE)

# explicitly set our default assay
# we want to use the original counts and not the integrated data
# DefaultAssay(integ10x) <- "RNA"

cluster0_conserved_markers <- FindConservedMarkers(integ10x,
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
  FindConservedMarkers(integ10x,
                       ident.1 = cluster,
                       grouping.var = "ITgeno",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    left_join(y = unique(annotations[, c("gene_name", "description")]),
              by = c("gene" = "gene_name")) %>%
    cbind(cluster_id = cluster, .)
}

# Iterate function across desired clusters
conserved_markers <- map_dfr(c(0:17), get_conserved)


# Extract top 10 markers per cluster
top10 <- conserved_markers %>% 
  mutate(avg_fc = (II_avg_log2FC + IT_avg_log2FC) /2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

# Visualize top 10 markers per cluster
View(top10)

# Plot interesting marker gene expression
FeaturePlot(object = integ10x, 
            features = c("FCGR2B"),
            order = TRUE,
            min.cutoff = 'q10', 
            split.by = 'ITgeno')

# Vln plot - cluster 
VlnPlot(object = integ10x, 
        features = c("FCGR2B"))

integ10x <- SetIdent(integ10x, value = 'ITgeno')

# save labeling seurat
write_rds(integ10x, 'data/cell/label10x.rds')

# subset tumor infiltrates
tumor10x <- subset(integ10x, Tissue == 'T')
write_rds(tumor10x, 'data/cell/tumor10x.rds')

# save top10 markers
write_csv(top10, 'results/cluster/10x-clustersMarker.csv')

# Determine differentiating markers for CD4+ T cell
# cd4_tcells <- FindMarkers(integ10x,
#                           ident.1 = 2,
#                           ident.2 = c(0,4,10))                  
# 
# # Add gene symbols to the DE table
# cd4_tcells <- cd4_tcells %>%
#   rownames_to_column(var = "gene") %>%
#   left_join(y = unique(annotations[, c("gene_name", "description")]),
#             by = c("gene" = "gene_name"))
# 
# # Reorder columns and sort by padj      
# cd4_tcells <- cd4_tcells[, c(1, 3:5,2,6:7)]
# 
# cd4_tcells <- cd4_tcells %>%
#   dplyr::arrange(p_val_adj) 
# 
# # View data
# View(cd4_tcells)
# 
# # show difference between our clusters and Zhang's clusters
# label10x@meta.data %>% 
#   ggplot(aes(x=Sub_Cluster, fill=integrated_snn_res.1)) + 
#   geom_bar() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   theme(plot.title = element_text(hjust=0.5, face="bold")) +
#   ggtitle("NCells")

