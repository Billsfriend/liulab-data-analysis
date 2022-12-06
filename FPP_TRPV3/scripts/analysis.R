library(Seurat)
library(tidyverse)
library(tidyseurat)
library(ggpubr)

read_rds('FPP_TRPV3/data/Dunlap-Rao_2022_CLE/sobj.rds') -> sobj

# cell type fraction
sobj@meta.data -> meta

sobj %>%
  unite('fine_label',
        c('singler_labels','seurat_clusters'),
        remove = FALSE) ->
  sobj

SetIdent(sobj, value = 'fine_label') -> sobj

sobj %>%
  filter(singler_labels == 'Keratinocytes')%>%
  ggplot(aes(orig.ident, fill = fine_label))+
  geom_bar()+
  coord_flip()+
  scale_fill_viridis_d(begin = 0.1, option = 'turbo')+
  labs_pubr()+
  labs(title = 'Heterogenity of keratinocytes in CLE-lesion/HC skin',
       x = 'group')

# cytokine in interest -------------

jzl_list <- read.table('CRC-I/ref/JZL_list')

sobj@assays$RNA@counts %>% rownames() -> genes

jzl_list %>%
  filter(V1 %in% genes) %>%
  pull(V1) -> cytokines_interest

find_interest_markers <- function(cell_type){
  try(FindMarkers(sobj,
                  subset.ident = cell_type,
                  group.by = 'orig.ident',
                  features = cytokines_interest,
                  ident.1 = 'Lupus_lesion') %>%
        add_column(cell_type = cell_type) %>%
        rownames_to_column('gene'))
}

type_name_list <- unique(sobj$fine_label)

map(type_name_list, find_interest_markers) %>%
  purrr::reduce(bind_rows) ->
  interest_deg_by_cell_type

interest_deg_by_cell_type %>%
  ggplot() +
  geom_point(aes(x = gene, y = cell_type, color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  labs_pubr()+
  ggtitle('DEG cytokines in interest in CLE-lesion skin vs HC')

# FPP pathway ----------
read.table('BD-MassSpec/data/fpp_genes.txt') -> fpp_pathway

upstream <- fpp_pathway$V1[1:11]
downstream <- fpp_pathway$V1[12:20]

cytokines_interest <- upstream

map(type_name_list, find_interest_markers) %>%
  purrr::reduce(bind_rows) ->
  interest_deg_by_cell_type

interest_deg_by_cell_type %>%
  ggplot() +
  geom_point(aes(x = gene, y = cell_type, color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  labs_pubr()+
  ggtitle('DEG upstream FPP synthesis in CLE-lesion skin vs HC')

VlnPlot(sobj, features = c('FDPS','PMVK'), idents = 'Keratinocytes_4', group.by  = 'orig.ident', pt.size = 0)

cytokines_interest <- downstream

map(type_name_list, find_interest_markers) %>%
  purrr::reduce(bind_rows) ->
  interest_deg_by_cell_type

interest_deg_by_cell_type %>%
  ggplot() +
  geom_point(aes(x = gene, y = cell_type, color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  labs_pubr()+
  ggtitle('DEG downstream FPP synthesis in CLE-lesion skin vs HC')


# inflammatory genes --------------
inflam_list <- c('IL6','CCL20','TSLP','CXCL8')

cytokines_interest <- inflam_list

map(type_name_list, find_interest_markers) %>%
  purrr::reduce(bind_rows) ->
  interest_deg_by_cell_type

interest_deg_by_cell_type %>%
  ggplot() +
  geom_point(aes(x = gene, y = cell_type, color = avg_log2FC, size = -log10(p_val_adj))) +
  scale_color_gradient2(low = 'blue', high = 'red') +
  labs_pubr()+
  ggtitle('DEG associated with inflammation in interest in CLE-lesion skin vs HC')

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
AddModuleScore(sobj, features = list(upstream, downstream), name = 'updown_fpp') -> sobj

VlnPlot(sobj, 'updown_fpp1', pt.size = 0, split.by = 'orig.ident') +
  ggtitle('FPP upstream pathway score')

VlnPlot(sobj, 'updown_fpp2', pt.size = 0, split.by = 'orig.ident') +
  ggtitle('FPP downstream pathway score')

Idents(sobj) <- 'singler_labels'

FindMarkers(sobj,
            subset.ident = 'Keratinocytes',
            group.by = 'fine_label',
            ident.1 = 'Keratinocytes_4') -> k4marker

k4marker %>%
  rownames_to_column('gene') %>%
  filter(p_val_adj < 0.05) %>%
  write_csv('FPP_TRPV3/results/k4marker.csv')

read_csv('FPP_TRPV3/results/k4marker.csv') -> k4marker

filter(k4marker, p_val_adj < 0.05) %>%
  slice_max(avg_log2FC, n = 4) %>%
  pull(gene) ->
  k4top_marker

VlnPlot(sobj,
        idents = 'Keratinocytes',
        group.by = 'fine_label',
        k4top_marker,
        pt.size = 0)

filter(k4marker, avg_log2FC >= 2) %>%
  pull(gene) %>%
  clusterProfiler::enrichGO(
    OrgDb = 'org.Hs.eg.db',
    keyType = 'SYMBOL',
    ont = 'BP'
  ) ->
  res_go_bp

res_go_bp@result -> t
