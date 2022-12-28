{
  library(Seurat)
  library(tidyverse)
  library(ggpubr)
  library(pheatmap)
  options(stringsAsFactors = F)
}

# load data
bgi_sobj <- readRDS("covid19/data/Final_nCoV_0716_upload_BGI.rds")

blish_sobj <- readRDS("covid19/data/blish_covid_NM2020.rds")

DefaultAssay(bgi_sobj) <- "RNA" # set default assay
DefaultAssay(blish_sobj) <- 'RNA'

blish_sobj <- subset(blish_sobj, Status == "Healthy") # only want the control group in this dataset
bgi_sobj <- subset(bgi_sobj, Stage == "Ctrl")

# use singler + monaco to re-annotate cell type
monaco <- celldex::MonacoImmuneData()

blish_sobj %>%
  as.SingleCellExperiment() %>%
  SingleR::SingleR(ref = monaco,
                   labels = monaco$label.fine,
                   clusters = blish_sobj$cell.type) ->
  blish_type_res

tibble(original = blish_type_res@rownames,
       singleR = blish_type_res$pruned.labels)->
  blishtype

bgi_sobj %>%
  as.SingleCellExperiment() %>%
  SingleR::SingleR(ref = monaco,
                   labels = monaco$label.fine,
                   clusters = bgi_sobj$cell_type) ->
  bgi_type_res


tibble(original = bgi_type_res@rownames,
       singleR = bgi_type_res$pruned.labels)->
  bgitype

blish_sobj@meta.data %>%
  left_join(blishtype, by = c('cell.type' = 'original')) ->
  blish_sobj@meta.data

bgi_sobj@meta.data %>%
  left_join(bgitype, by = c('cell_type' = 'original')) ->
  bgi_sobj@meta.data

Idents(bgi_sobj) <- 'singleR'
Idents(blish_sobj) <- 'singleR'

table(Idents(bgi_sobj))

# T & B cell marker -------------
t_features <- c('CD3D','CD3E','CD3G','CD247','TRAC','TRBC1','TRBC2','TRDC','TRGC1','TRGC2')

b_features <- c('CD79A','CD79B','IGKC','IGLC2','IGLC3','IGLC6','IGLC7','IGHA1','IGHA2','IGHG1','IGHG2','IGHG3','IGHG4','IGHD','IGHE','IGHM')

cross_tibble <- expand_grid(
  sobj = c('blish', 'bgi'),
  feature = c('t_features', 'b_features'))

find_hook <- function(hook){
  switch (hook,
          'bgi' = bgi_sobj,
          'blish' = blish_sobj,
          'b_features' = b_features,
          't_features' = t_features
  )
}

# violin plot-----------
plot_hooked_violin <- function(sobj,features){
  VlnPlot(find_hook(sobj),
          features = find_hook(features),
          pt.size = 0) &
    theme(axis.title.x = element_blank(),
          text = element_text(size = 5),
          axis.text = element_text(size = 5))}

map2(cross_tibble$sobj,
     cross_tibble$feature,
     plot_hooked_violin) ->
  violinp

violinp[[1]]

# dot plot --------------
plot_hooked_dotplot <- function(sobj, features){
  DotPlot(find_hook(sobj),
          features = find_hook(features),
          scale = FALSE,
          cols = c('blue','red'),
          cluster.idents = TRUE) +
    RotatedAxis()+
    ylab('Cell types')}

map2(cross_tibble$sobj,
     cross_tibble$feature,
     plot_hooked_dotplot) ->
  dotplotp

dotplotp[[1]]

ggsave("DE_cells/figures/bgi_b_receptors_dotplot.pdf",
       width = 20,
       height = 15,
       units = 'cm')

# heatmap ------------
percent_above0 <- function(x){
  PercentAbove(x, threshold = 0)
}

plot_hooked_heatmap <- function(sobj, features){
  pos_percent_for_types <- function(type_name){
    find_hook(sobj) %>%
      subset(ident = type_name) %>%
      FetchData(vars = find_hook(features)) %>%
      lapply(percent_above0) %>%
      as.data.frame() %>%
      add_column(cell_type = type_name)
  }
  
  find_hook(sobj) %>%
    Idents() %>%
    unique() %>%
    lapply(pos_percent_for_types) %>%
    purrr::reduce(bind_rows) ->
    heat_mat
  
  heat_mat %>%
    column_to_rownames('cell_type') %>%
    mutate(across(everything(), ~formatC(.*100, format = 'f', digits = 2))) %>%
    mutate(across(everything(), ~paste0(.,'%'))) ->
    percent_mat_chr
  
  heat_mat %>%
    column_to_rownames('cell_type') %>%
    t() %>%
    pheatmap(
      cluster_row = FALSE,
      cluster_cols = FALSE,
      display_numbers = t(percent_mat_chr),
      number_color = "black",
      legend = FALSE,
      fontsize = 6,
      fontsize_number = 7
    ) %>%
    ggplotify::as.ggplot()}

map2(cross_tibble$sobj,
     cross_tibble$feature,
     plot_hooked_heatmap) ->
  heatmapp

add_column(cross_tibble,
           heatmap = heatmapp,
           dotplot = dotplotp,
           violin = lapply(violinp, ggplotify::as.ggplot)) ->
  fig_tibble

fig_tibble$heatmap[[3]] +
  fig_tibble$violin[[3]] +
  fig_tibble$heatmap[[1]] +
  fig_tibble$violin[[1]] +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18))

ggsave("DE_cells/figures/ms_t_marker_fig1.pdf",
       width = 36,
       height = 30,
       units = 'cm')

fig_tibble$heatmap[[4]] +
  fig_tibble$violin[[4]] +
  fig_tibble$heatmap[[2]] +
  fig_tibble$violin[[2]] +
  patchwork::plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(size = 18))

ggsave("DE_cells/figures/ms_b_marker_fig2.pdf",
       width = 36,
       height = 30,
       units = 'cm')


