library(data.table)
library(Seurat)
library(harmony)
library(pheatmap)
library(enrichR)
library(SingleR)
library(BiocParallel)
library(tidyseurat)
library(tidyverse)

# read raw data -------------
count_mtx <- fread('data/GSE162183_psoriasis_gene_counts_matrix.tab.gz')

count_mtx[1:5, 1:5]

column_to_rownames(count_mtx, 'V1') -> count_mtx

CreateSeuratObject(count_mtx) -> sobj

sobj$any <- 1

VlnPlot(sobj, c('nFeature_RNA', 'nCount_RNA'), group.by = 'any')

sobj@meta.data -> metadata

metadata %>%
  mutate(group = case_when(
    str_detect(orig.ident, 'Ctrl') ~ 'Ctrl',
    TRUE ~ 'Psor')) -> 
  sobj@meta.data

SetIdent(sobj, value = 'group') -> sobj

# pre-process -------------
sobj %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA() -> sobj

DimPlot(sobj)

RunHarmony(sobj, group.by.vars = 'orig.ident') -> sobj

sobj %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity() ->
  sobj

DimPlot(sobj)

# cluster identify ---------
human_cell_ref <- read_rds('CRC-I/ref/HumanPrimaryCellAtlas.rds')

per_clust_idented <- sobj %>%
  as.SingleCellExperiment() %>%
  SingleR(
    ref = human_cell_ref,
    labels = human_cell_ref$label.main,
    clusters = sobj$seurat_clusters,
    de.method = "wilcox",
    BPPARAM = MulticoreParam()
  )

table(per_clust_idented$pruned.labels)

plotScoreHeatmap(per_clust_idented,
                 show_colnames = TRUE,
                 show.pruned = TRUE)

match(sobj@meta.data$seurat_clusters,
      rownames(per_clust_idented)) -> match_index

sobj$singler_labels <- per_clust_idented$pruned.labels[match_index]

sobj %>%
  mutate(labels = case_when(
    singler_labels == 'iPS_cells' ~ 'Fibroblasts',
    singler_labels == 'Astrocyte' ~ 'Neurons',
    singler_labels == 'BM & Prog.' ~ 'NK cells',
    TRUE ~ singler_labels
        )) -> sobj

SetIdent(sobj, value = 'labels') -> sobj

# save and read -----------
write_rds(sobj, 'data/other/psoriatis0909.rds')

sobj <- read_rds('data/other/psoriatis0909.rds')

# DEG by cell type ---------
unique(sobj$labels) -> type_name_list

jzl_list <- read.table('CRC-I/ref/JZL_list')

find_interest_markers <- function(cell_type){
  try(FindMarkers(sobj,
                  #subset.ident = cell_type,
                  group.by = 'group',
                  features = jzl_list$V1,
                  ident.1 = 'Psor') %>%
        add_column(cell_type = cell_type) %>%
        rownames_to_column('gene'))
}

map(type_name_list, find_interest_markers) %>%
  reduce(bind_rows) ->
  interest_deg_by_cell_type

interest_deg_by_cell_type %>%
  filter(gene %in% jzl_list$V1) %>%
  ggplot() +
  geom_point(aes(x = gene, y = cell_type, color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  ggtitle('cytokines in interest')

# visualization mission -----------
read.table('BD-MassSpec/data/fpp_genes.txt') -> fpp_pathway

VlnPlot(sobj, 'fpp_program1', pt.size = 0) +
  ggtitle('FPP pathway score')

DimPlot(sobj)

# upstream/downstream FPP
upstream <- fpp_pathway$V1[1:11]
downstream <- fpp_pathway$V1[12:20]

filter(sobj, labels == 'Keratinocytes') %>%
DotPlot(features = upstream, group.by = 'orig.ident') +
  labs(y = 'group',
       title = 'pathway upstream FPP')+
  coord_flip() -> p1

filter(sobj, labels == 'Keratinocytes') %>%
DotPlot(features = downstream, group.by = 'orig.ident') +
  labs(y = 'group',
       title = 'pathway downstream FPP')+
  coord_flip() -> p2

p1 + p2


# inflammatory genes
inflam_list <- c('IL6','CCL20','TSLP','CXCL8')
DotPlot(sobj, features = inflam_list, group.by = 'orig.ident')+
  labs(y = 'group',                                    title = 'cytokine expression')+
  coord_flip()

VlnPlot(sobj, features = downstream, pt.size = 0)

AverageExpression(sobj, features = inflam_list, group.by = 'orig.ident') -> inflam_tbl

pheatmap(inflam_tbl$RNA,
         cluster_cols = FALSE,
         angle_col = 0)

inflam_tbl$RNA %>%
  as_tibble(rownames = 'gene') %>%
  pivot_longer(2:7, names_to = 'sample', values_to = 'mean_expression') %>%
  mutate(group = str_sub(sample, 1, 4))->
  mean_inflam_expr

mean_inflam_expr %>%
  filter(group == 'Ctrl' & gene == 'CCL20')%>%
  pull(mean_expression) -> ctrl_ccl20

mean_inflam_expr %>%
  filter(group == 'Psor' & gene == 'CCL20')%>%
  pull(mean_expression) -> psor_ccl20

t.test(ctrl_ccl20, psor_ccl20)

ggplot(mean_inflam_expr, aes(x = group, y = mean_expression)) +
  geom_jitter(aes(shape = group), width = 0.03) +
  stat_summary(geom = "crossbar", fun = mean, colour = "red", width = 0.2)+
  facet_wrap(vars(gene), scales = 'free')

# use this integrated function to do quick pathway enrichment!
enrichR::listEnrichrSites()
Seurat::DEenrichRPlot(sobj,
                      ident.1 = 'Psor',
                      enrich.database = 'GO_Biological_Process_2021',
                      max.genes = 1000)

# calculate program score like papers --------
Seurat::AddModuleScore(sobj, features = fpp_pathway, name = 'fpp_program') -> sobj




