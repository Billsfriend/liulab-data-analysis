{
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(cowplot)
  library(pheatmap)
  options(stringsAsFactors = F)
}

# load data
Human_PBMC <- readRDS("covid19/data/Final_nCoV_0716_upload_BGI.RDS")
DefaultAssay(Human_PBMC) <- "RNA" # set default assay
pbmc <- subset(Human_PBMC, idents = "Ctrl") # only want the control group in this dataset

my_level <-
  c(
    "Naive B cells",
    "Memory B cells",
    "Plasma",
    "Cycling Plasma",
    "Naive T cells",
    "Activated CD4 T cells",
    "Cytotoxic CD8 T cells",
    "Cycling T cells",
    "MAIT",
    "NKs",
    "XCL+ NKs",
    "Monocytes",
    "Megakaryocytes",
    "DCs",
    "Stem cells"
  )

table(pbmc$cell_type) # make it clear to watch

Idents(pbmc) <- factor(pbmc$cell_type, levels = my_level)

table(Idents(pbmc)) # make sure only get the control group

#T cell marker -------------
features <- c('CD3D','CD3E','CD3G','CD247','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2')

# violin plot-----------
VlnPlot(pbmc, 
        features = features, 
        pt.size = 0,
        group.by = 'cell_type') &
  theme(axis.title.x = element_blank(),
        text = element_text(size = 6),
        axis.text = element_text(size = 6)) -> p1

ggsave("DE_cells/figures/bgi-t_receptors_violin.pdf",
       width = 18,
       height = 18,
       units = 'cm')


# dot plot --------------
DotPlot(pbmc, features = features, group.by = 'cell_type') +
  RotatedAxis() +
  scale_color_gradient(low = 'blue', high = 'red')
ggsave("DE_cells/figures/bgi-t_receptors_dotplot.pdf",
       width = 12,
       height = 10,
       units = 'cm')


# heatmap ------------
percent_above0 <- function(x){
  PercentAbove(x, threshold = 0)
  }

pos_percent_for_types <- function(type_name){
  pbmc %>%
    subset(cell_type == type_name) %>%
    FetchData(vars = features) %>%
    lapply(percent_above0) %>%
    as.data.frame() %>%
    add_column(cell_type = type_name)
}

lapply(my_level, pos_percent_for_types) %>%
  purrr::reduce(bind_rows) ->
  heat_mat

heat_mat %>%
  column_to_rownames('cell_type') %>%
  mutate(across(everything(), ~.*100)) %>%
  mutate(across(everything(), as.character)) %>%
  mutate(across(everything(), ~str_trunc(., width = 5, side = 'right', ellipsis = '%'))) ->
  percent_mat_chr

heat_mat %>%
  column_to_rownames('cell_type') %>%
  pheatmap(
    cluster_row = FALSE, 
    cluster_cols = FALSE,
    display_numbers = percent_mat_chr,
    number_color = "black",
    legend = FALSE
) -> p2

ggplotify::as.ggplot(p2) + p3 +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18))

ggsave("DE_cells/figures/ms_t_marker_fig1.pdf", 
       width = 36, 
       height = 22,
       units = 'cm')


