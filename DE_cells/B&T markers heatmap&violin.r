{
  library(Seurat)
  library(tidyverse)
  library(ggpubr)
  library(clusterProfiler)
  library(pheatmap)
  options(stringsAsFactors = F)
}

# load data
Human_sobj <- readRDS("covid19/data/Final_nCoV_0716_upload_BGI.rds")

Human_sobj <- readRDS("covid19/data/blish_covid_NM2020.rds")

DefaultAssay(Human_sobj) <- "RNA" # set default assay

sobj <- subset(Human_sobj, Status == "Healthy") # only want the control group in this dataset

Idents(sobj) <- 'cell.type'

table(sobj$cell_type) # make it clear to watch

table(Idents(sobj)) # make sure only get the control group

my_level <- unique(sobj$cell.type)

#T cell marker -------------
features <- c('CD3D','CD3E','CD3G','CD247','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2')

# B cell markers --------------
features <- c('CD79A','CD79B','IGKC','IGLC2','IGLC3','IGLC6','IGLC7','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM')


# violin plot-----------
VlnPlot(sobj, 
        features = features, 
        pt.size = 0,
        group.by = 'cell.type',
        ncol = 3) &
  theme(axis.title.x = element_blank(),
        text = element_text(size = 5),
        axis.text = element_text(size = 5)) -> p1

ggsave("DE_cells/figures/bgi-t_receptors_violin.pdf",
       width = 18,
       height = 18,
       units = 'cm')


# dot plot --------------
DotPlot(sobj, features = features) +
  RotatedAxis() +
  scale_color_gradient(low = 'blue', high = 'red')
ggsave("DE_cells/figures/blish-t_receptors_dotplot.pdf",
       width = 12,
       height = 10,
       units = 'cm')

# heatmap ------------
percent_above0 <- function(x){
  PercentAbove(x, threshold = 0)
  }

pos_percent_for_types <- function(type_name){
  sobj %>%
    subset(ident = type_name) %>%
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
  mutate(across(everything(), ~formatC(.*100, format = 'f', digits = 2))) %>%
  mutate(across(everything(), ~paste0(.,'%'))) ->
  percent_mat_chr

heat_mat %>%
  column_to_rownames('cell_type') %>%
  pheatmap(
    cluster_row = FALSE, 
    cluster_cols = FALSE,
    display_numbers = percent_mat_chr,
    number_color = "black",
    legend = FALSE,
    fontsize = 6,
    fontsize_number = 6
) -> p2

ggplotify::as.ggplot(p2) + ggplotify::as.ggplot(p1) +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18))

ggsave("DE_cells/figures/ms_t_marker_fig1.pdf", 
       width = 30, 
       height = 12,
       units = 'cm')


