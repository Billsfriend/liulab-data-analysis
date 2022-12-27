library(tidyverse)
library(ggpubr)

read_csv('covid19/data/wu-chen/Protein_matrix_tidy.csv') -> prot_tidy

sample_info <- read_tsv('covid19/data/wu-chen-samples')

prot_tidy %>%
  filter(IGHG1_G396R != 'GR') %>%
  pivot_longer(where(is.numeric), names_to = 'gene', values_to = 'value') %>%
  group_by(IGHG1_GR, gene) %>%
  select(gene, IGHG1_GR, value) %>%
  summarise(value = list(value)) %>%
  pivot_wider(names_from = IGHG1_GR, values_from = value) %>%
  rowwise() %>%
  mutate(log2fc = log2(
    mean(unlist(RR))/mean(unlist(GG_GR))
  )) %>%
  mutate(p.value = wilcox.test(
    unlist(GG_GR), unlist(RR)
  )$p.value) %>%
  select(gene, p.value, log2fc) -> res

res %>%
  ungroup() %>%
  mutate(p.adj = p.adjust(p.value, 'BH')) %>%
  filter(p.value < 0.05) -> sig_res

write_csv(sig_res, 'covid19/results/wu-chen-prot-sig.csv')

# go analysis -----------
sig_res %>%
  filter(log2fc > 0) %>%
  pull(gene) %>%
  clusterProfiler::enrichGO(
    OrgDb = 'org.Hs.eg.db',
    keyType = 'SYMBOL',
    ont = 'BP'
  ) -> res_go

res_go@result %>%
  as_tibble() %>%
  filter(p.adjust<0.05) -> res_go_enrich

write_csv(res_go_enrich, 'covid19/results/wu-chen-prot-go-RRvsGG.csv')

save_pdf <- function(filename){
  ggsave(filename = filename,
         width = 24,
         height = 14,
         units = 'cm')
}

res %>%
  mutate(type = case_when(
    p.value < 0.05 & log2fc > 0 ~ 'Up in RR',
    p.value < 0.05 & log2fc < 0 ~ 'Down in RR',
    TRUE ~ 'NS'
  )) %>%
  ggplot(aes(log2fc, -log10(p.value), color = type))+
  geom_point() +
  geom_hline(yintercept = 1.3, linetype = 'dashed')

read_delim('covid19/results/wu-chen-prot-RRvsGG-digest.txt') -> prot_go

prot_go %>%
  mutate(Description = str_wrap(Description, 40)) %>%
  mutate(Description = reorder(Description, Count)) %>%
  ggplot(aes(Description, GeneRatio, size = Count, color = p.adjust)) +
  geom_point()+
  expand_limits(y = c(0,0.28))+
  scale_size(breaks = 2:4)+
  theme_bw()+
  labs_pubr()+
  scale_color_gradient(low = 'red', high = 'blue')+
  coord_flip()+
  labs(x = NULL,
       y = 'Gene ratio',
       title = 'Enriched pathways in\n PBMC proteome',
       subtitle = 'RR vs GG')+
  guides(colour = guide_colorbar(order = 2),
         size = guide_legend(order = 1))
  
ggsave('covid19/figures/prot-go.pdf',
       width = 16,
       height = 11,
       units = 'cm')

# boxplot of prot abundance -----------
prot_tidy %>%
  pivot_longer(1:634, names_to = 'gene', values_to = 'value') %>%
  filter(gene %in% c('C1QA','F2','HRG') & !is.na(IGHG1_G396R)) %>%
  ggplot(aes(IGHG1_G396R, value, color = IGHG1_G396R)) +
  geom_boxplot() +
  geom_jitter(width = 0.1, height = 0.1) +
  facet_wrap(~gene, scales = 'free')+
  labs(y = 'log2(abundance+1)', title = 'GO pathway: Humoral immune response') +
  theme_pubr() +
  labs_pubr()+
  stat_compare_means(comparisons = list(c('GR','RR'),c('GG','RR')), label = 'p.signif', step.increase = 0.3)+
  scale_color_manual(values = c('black','blue','red'))+
  scale_y_continuous(expand = expansion(mult = c(0.1,0.2)))

save_pdf('covid19/figures/prot-boxplot-humoral immune response.pdf')

prot_tidy %>%
  pivot_longer(1:634, names_to = 'gene', values_to = 'value') %>%
  filter(gene %in% c('C1QA','F2','HRG') & !is.na(IGHG1_G396R)) %>%
  select(id1, IGHG1_G396R, gene, value) %>%
  pivot_wider(names_from = 'gene', values_from = 'value') -> xls_prot

write_csv(xls, 'covid19/results/source_data_prot_go_Humoral_response.csv')

# GO GSEA ------------
res$log2fc -> gse_list

names(gse_list) <- res$gene

gse_list <- sort(gse_list, decreasing = TRUE)

clusterProfiler::gseGO(gse_list,
                       OrgDb = 'org.Hs.eg.db',
                       keyType = 'SYMBOL') ->
  res_gse

# concordance of mRNA and protein ------------
prot_sig <- read_csv('covid19/results/wu-chen-prot-sig.csv')

rna_sig <- read_delim('covid19/results/wu-chen-sdeg-RRvsGG_GR.tsv')

prot_sig %>%
  filter(log2fc > 0) %>%
  pull(gene) -> up_prot

rna_sig %>%
  filter(logFC > 0 & SYMBOL %in% up_prot) ->
  up_both

clusterProfiler::enrichGO(up_both$SYMBOL,
                          OrgDb = 'org.Hs.eg.db',
                          keyType = 'SYMBOL',
                          ont = 'BP') ->
  res_go_both

res_go_both@result %>%
  as_tibble() %>%
  filter(p.adjust < 0.05) ->
  topGo_both
  
write_csv(topGo_both, 'covid19/results/wu-chen-2omics-go.csv')

# manually selected pathways ---------
read_delim('covid19/results/wu-chen-2omics-digest.txt') %>%
  select(c(Description, geneID))%>%
  separate_rows(geneID) %>%
  distinct(geneID, .keep_all = TRUE) -> digest_2omics

cpm <- read_delim('covid19/results/wu-chen-cpm.tsv')

cpm %>%
  filter(SYMBOL %in% digest_2omics$geneID) %>%
  column_to_rownames('SYMBOL') %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('run_accession') %>%
  left_join(sample_info) %>%
  pivot_longer(cols = where(is.numeric), names_to = 'gene', values_to = 'cpm') -> digest_tidy

digest_tidy %>%
  group_by(gene, IGHG1_G396R) %>%
  summarize(median = median(cpm)) %>%
  pull(IGHG1_G396R) -> genotype_list

digest_tidy %>%
  group_by(gene, IGHG1_G396R) %>%
  summarize(median = median(cpm)) %>%
  group_by(gene) %>%
  summarize(zscore = (median - mean(median))/sd(median)) %>%
  add_column(IGHG1_G396R = genotype_list, type = 'mRNA') ->
  res_2omics

prot_tidy %>%
  filter(!is.na(IGHG1_G396R)) %>%
  select(one_of(digest_2omics$geneID),IGHG1_G396R) %>%
  pivot_longer(cols = where(is_numeric), names_to = 'gene',
               values_to = 'prot') %>%
  group_by(IGHG1_G396R, gene) %>%
  summarize(median = median(prot)) %>%
  group_by(gene) %>%
  summarize(zscore = (median - mean(median))/sd(median)) %>%
  add_column(IGHG1_G396R = genotype_list, type = 'protein') %>%
  bind_rows(res_2omics)->
  final_2omics

final_2omics %>%
  mutate(pathway = case_when(
    gene %in% c('F2','LTF') ~ 'antimicrobial humoral response',
    TRUE ~ 'response to interferon-gamma'
  )) %>%
  ggplot(aes(IGHG1_G396R, zscore, color = type, group = type))+
  geom_point()+
  geom_path()+
  facet_wrap(~pathway+gene, scales = 'free', labeller = 'label_both')+
  theme_pubr()+
  labs_pubr()+
  labs(y = 'scaled median expression',
       title = 'Concordance between upregulated mRNA and protein in RR genotype PBMC')

ggsave('covid19/figures/wu-chen-2omics-up-in-RR.pdf',
       width = 20,
       height = 15,
       units = 'cm')  
