library(edgeR)
library(tidyverse)
library(ggpubr)

data <- read_tsv('covid19/data/wu-chen-counts.tsv')

# process metadata -------------
sample_info <- read_tsv('covid19/data/wu-chen-samples') 

my_pal <- c('black', 'blue', 'red')

save_pdf <- function(filename, width = 24, height = 14){
  ggsave(filename = filename,
         width = width,
         height = height,
         units = 'cm')
}

sample_info %>%
  ggplot(aes(clinic, fill = IGHG1_G396R)) +
  geom_bar() +
  theme_pubr() +
  labs_pubr()+
  labs(x = 'disease stage') +
  scale_fill_manual(values = my_pal)+
  coord_flip() -> p1
  
sample_info %>%
  ggplot(aes(clinic, fill = IGHG1_G396R)) +
  geom_bar(position = 'fill') +
  theme_pubr() +
  labs_pubr()+
  labs(x = 'disease stage', y = 'proportion') +
  scale_fill_manual(values = my_pal)+
  coord_flip()-> p2

p1 / p2

save_pdf('covid19/figures/GR_severity_test.pdf', width = 16, height = 9)

data %>%
  column_to_rownames('SYMBOL') ->
  data

sample_info %>%
  filter() %>%
  pull(run_accession) ->
  sample_list

# generate DGEList object -------------
data %>%
  select(one_of(sample_list)) %>%
  DGEList() -> edglist

edglist$samples %>%
  rownames_to_column('run_accession') %>%
  left_join(sample_info) %>%
  column_to_rownames('run_accession') -> sample_data

sample_data$IGHG1_GR %>%
  as_factor() %>%
  relevel(ref = 'GG_GR') -> sample_data$group

# sample_data$clinic %>%
#   as_factor() -> sample_data$clinic

edglist$samples <- sample_data

#design <- model.matrix(~sample_data$clinic + sample_data$group)

#rownames(design) <- rownames(sample_data)

keep <- filterByExpr(edglist)
edglist <- edglist[keep, , keep.lib.sizes = FALSE]

edglist %>%
  calcNormFactors() %>%
  estimateDisp() -> edglist

# plotMDS(edglist, labels = edglist$samples$group)

# for more aggressive search for DEG:
edglist %>%
  glmFit() %>%
  glmLRT() -> res_glmLRT

res_glmLRT %>%
  decideTests() %>%
  summary()

# for more safe and conservative
edglist %>%
  glmQLFit() %>%
  glmQLFTest() -> res_glmQL

res_glmQL %>%
  decideTests() %>%
  summary()

# alternative significance threshold
edglist %>%
  glmFit() %>%
  glmTreat() -> res_glmTreat

res_glmTreat %>%
  decideTests() %>%
  summary()

edglist %>%
  glmQLFit() %>%
  glmTreat() -> res_glmQLTreat

res_glmQLTreat %>%
  decideTests() %>%
  summary()

# visualization ----------
res_glmLRT$table %>%
  rownames_to_column('SYMBOL') -> sdeg

sdeg$padj <- p.adjust(sdeg$PValue, method = 'BH')

sdeg %>%
  mutate(type = case_when(
    padj < 0.05 & logFC > 1 ~ 'Up in RR',
    padj < 0.05 & logFC < -1 ~ 'Down in RR',
    TRUE ~ 'NS'
  )) -> sdeg

EnhancedVolcano::EnhancedVolcano(
  sdeg,
  x = 'logFC',
  y = 'PValue',
  lab = sdeg$SYMBOL,
  title = NULL,
  subtitle = NULL,
  caption = NULL)

write_tsv(sdeg, 'covid19/results/wu-chen-sdeg-RRvsGG_GR.tsv')

# compute cpm --------------
cpm(edglist) %>%
  as.data.frame() %>%
  rownames_to_column('ENTREZID') %>%
  right_join(sdeg[,1:2]) %>%
  relocate(SYMBOL) %>%
  filter(!is.na(SYMBOL)) %>%
  select(-ENTREZID) -> cpm

write_tsv(cpm, 'covid19/results/wu-chen-cpm.tsv')

# GO analysis ------------
if (!exists('sdeg')){
  sdeg <- read_tsv('covid19/results/wu-chen-sdeg.tsv')
}

sdeg %>%
  filter(padj < 0.05 & logFC > 0) %>%
  pull(SYMBOL) %>%
  clusterProfiler::enrichGO(
    OrgDb = 'org.Hs.eg.db',
    ont = 'BP',
    keyType = 'SYMBOL'
  ) -> res_go

res_go@result %>%
  as_tibble() %>%
  filter(p.adjust < 0.05) %>%
  select(-c(BgRatio, pvalue, qvalue)) %>%
  mutate(GeneRatio = str_extract(GeneRatio, '[:digit:]+$') %>% as.numeric()) %>%
  mutate(GeneRatio = Count/GeneRatio) -> res_topGo

write_csv(res_topGo, 'covid19/results/wu-chen-go-RRvsGG_GR.csv')

# GO GSEA -----------
sdeg$logFC -> logfc_set

names(logfc_set) <- sdeg$SYMBOL

sort(logfc_set, decreasing = TRUE) %>%
  clusterProfiler::gseGO(OrgDb = 'org.Hs.eg.db',
                       keyType = 'SYMBOL') ->
  res_gsea

res_gsea@result %>%
  as_tibble() %>%
  arrange(NES)-> res_topGse

res_topGse %>%
  slice_max(order_by = NES, n = 6) %>%
  bind_rows(slice_min(res_topGse, order_by = NES, n = 6)) %>%
  mutate(Description = reorder(Description, NES)) %>%
  ggplot(aes(Description, NES, fill = p.adjust))+
  geom_col()+
  coord_flip()+
  scale_fill_gradient(low = 'red', high = 'blue')+
  theme_bw()+
  labs_pubr()
  
write_excel_csv(res_topGse, 'covid19/results/wu-chen-gsea-RRvsGR.csv')

# after manually select pathway terms -----------
go_digest <- read_delim('covid19/results/wu-chen-go-RRvsGG-digest.txt')

go_digest %>%
  mutate(Description = str_wrap(Description, 40)) %>%
  mutate(Description = reorder(Description, GeneRatio))-> go_digest

ggplot(go_digest, aes(Description, GeneRatio, color = p.adjust, size = Count)) +
  geom_point()+
  coord_flip()+
  theme_bw()+
  labs_pubr()+
  scale_color_gradient(low = 'red', high = 'blue')+
  scale_size()+
  guides(colour = guide_colorbar(order = 2),
         size = guide_legend(order = 1))+
  labs(x = NULL,
       y = 'Gene ratio',
       title = 'Enriched pathways in\n PBMC transcriptome',
       subtitle = 'RR vs GG')

ggsave('covid19/figures/RNA-go-RRvsGG.pdf',
       width = 16,
       height = 11,
       units = 'cm')

# cpm boxplot of enriched go -----------
if (!exists('cpm')){
  cpm <- read_delim('covid19/results/wu-chen-cpm.tsv')
}

go_digest %>%
  select(c(Description, geneID)) %>%
  separate_rows(geneID) %>%
  left_join(sdeg, by = c('geneID' = 'SYMBOL')) %>%
  group_by(Description) %>%
  slice_min(order_by = padj, n = 9, with_ties = FALSE) -> digest_pathway
  
for (name in c('mucosal','heterop')) {
  digest_pathway %>%
  filter(str_detect(Description, 'B cell activa')) %>%
  slice_min(order_by = padj, n = 9) ->
  one_pathway

pathway_name <- one_pathway$Description[1]

cpm %>%
  filter(SYMBOL %in% one_pathway$geneID) %>%
  column_to_rownames('SYMBOL') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('run_accession') %>%
  left_join(sample_info) %>%
  pivot_longer(cols = where(is.numeric), names_to = 'gene', values_to = 'cpm') -> digest_tidy

digest_tidy %>%
  ggplot(aes(IGHG1_G396R, log(cpm+1), color = IGHG1_G396R))+
  geom_boxplot()+
  geom_jitter(width = 0.1, height = 0.1)+
  facet_wrap(~gene, scales = 'free')+
  theme_pubr()+
  labs_pubr()+
  scale_color_manual(values = my_pal) +
  ggtitle(str_glue('GO pathway: {pathway_name}')) +
  stat_compare_means(comparisons = list(c('GR','RR'),c('GG','RR')),label = 'p.signif', step.increase = 0.3)+
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))

ggsave(str_glue('covid19/figures/go-cpm-{pathway_name}.pdf'),
       width = 24,
       height = 18,
       units = 'cm')}

digest_tidy %>%
  select(IGHG1_G396R, id1, gene, cpm) %>%
  mutate(cpm = log(cpm+1)) %>%
  pivot_wider(names_from = gene, values_from = cpm) ->
  xls

write_csv(xls, 'covid19/results/source_data_go_B_activate.csv')

# cibersortx results ---------
read_tsv('covid19/results/CIBERSORTx_Job8_Adjusted.txt') -> nsclc

nsclc %>%
  rename(run_accession = Mixture) %>%
  left_join(sample_info) -> nsclc

pivot_longer(nsclc, cols = 2:7, names_to = 'cell_type', values_to = 'fraction') -> nsclc_tidy

# 3 group
ggplot(nsclc_tidy, aes(IGHG1_G396R, fraction, color = IGHG1_G396R)) +
  stat_summary(geom = 'col', fun = 'mean', fill = 'white', linewidth = 1) +
  stat_summary(geom = 'errorbar', fun.data = 'mean_cl_boot', width = 0.5)+
  geom_jitter(width = 0.1, height = 0)+
  facet_wrap(~cell_type, scales = 'free') +
  stat_compare_means(comparisons = list(c('GG','RR')),label.y.npc = 0.8, label = 'p.format')+
  theme_pubr()+
  labs_pubr()+
  scale_color_manual(values = my_pal)+
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))

ggsave('covid19/figures/cell_fraction.pdf',
       width = 24,
       height = 14,
       units = 'cm')

# 2 group
ggplot(nsclc_tidy, aes(IGHG1_GR, fraction, color = IGHG1_GR)) +
  stat_summary(geom = 'col', fun = 'mean', fill = 'white', linewidth = 1) +
  stat_summary(geom = 'errorbar', fun.data = 'mean_cl_boot', width = 0.5)+
  geom_jitter(width = 0.1, height = 0)+
  facet_wrap(~cell_type, scales = 'free') +
  stat_compare_means(comparisons = list(c('GG_GR','RR')),label.y.npc = 0.8, label = 'p.format')+
  theme_pubr()+
  labs_pubr()+
  scale_color_manual(values = c('skyblue','red'))+
  scale_y_continuous(expand = expansion(mult = c(0,0.2)))

ggsave('covid19/figures/cell_fraction_2.pdf',
       width = 24,
       height = 14,
       units = 'cm')
