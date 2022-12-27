library(DESeq2)
library(tidyverse)

data <- read_tsv('../covid19/data/wu-chen-counts.tsv')

sample_info <- read_tsv('../covid19/data/wu-chen-samples')

data %>%
  column_to_rownames('SYMBOL') -> data 

sample_info %>%
  filter() ->
  subset_info

data %>%
  select(one_of(subset_info$run_accession)) %>%
  DESeqDataSetFromMatrix(
  colData = subset_info,
  design = ~ clinic + IGHG1_GR) -> dds

# pre-filtering data of too low counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

# relevel to specify the control group (reference)
levels(dds$IGHG1_GR)
dds$IGHG1_GR <- relevel(dds$IGHG1_GR, ref = "GG_GR")

# rlog() may take a long time with 50 or more samples,
# vst() is a much faster transformation

# set parallel workers reduce time cost from 188s to 96s for 63*24200 matrix
dds <- DESeq(dds, parallel = TRUE)
res <- results(dds)
summary(res)

as.data.frame(res) %>%
  rownames_to_column('SYMBOL') -> deg_tibble

write_csv(deg_tibble, '../covid19/results/wu-chen-deseq2-ALL.csv')

read_csv('../covid19/results/wu-chen-deseq2-ALL.csv') -> deg_tibble

deg_tibble %>%
  mutate(type = case_when(
    log2FoldChange > 1 & pvalue < 0.01 ~ 'Enriched in RR',
    log2FoldChange < -1 & pvalue < 0.01 ~ 'Depleted in RR',
    TRUE ~ 'NS'
  )) %>%
  mutate(sig_symbol = case_when(
    SYMBOL %in% deg_b_acti_list$SYMBOL ~ 'red',
    TRUE ~ 'grey'
  )) -> deg_tibble

ggplot(deg_b_acti_list, aes(log2FoldChange, -log10(pvalue))) +
  geom_point() +
  geom_vline(xintercept = c(1, -1), linetype = 'dashed') +
  geom_hline(yintercept = 2, linetype = 'dashed') +
  #scale_color_manual(values = c('red', 'blue', 'grey')) +
  ggpubr::theme_pubr()

# clusterProfiler -----------
deg_tibble %>%
  filter(log2FoldChange > 0, pvalue < 0.01) %>%
  pull(SYMBOL) -> up_list

deg_tibble %>%
  filter(log2FoldChange < 0, pvalue < 0.01) %>%
  pull(SYMBOL) -> down_list

clusterProfiler::enrichGO(up_list,
         OrgDb = 'org.Hs.eg.db',
         keyType = 'SYMBOL',
         ont = 'BP',
         universe = deg_tibble$SYMBOL,
         readable = TRUE) -> ora_go

ora_go@result %>%
  filter(p.adjust < 0.1) %>%
  separate(GeneRatio, into = c('count','size')) %>%
  mutate(GeneRatio = Count / as.numeric(size)) ->
  ora_go_res

ora_go_res %>%
  ggplot(aes(Description, GeneRatio, color = p.adjust, size = Count))+
  geom_point() +
  coord_flip() +
  ggpubr::labs_pubr() +
  scale_color_viridis_c()

clusterProfiler::dotplot(ora_go)

cnetplot(ora_go)

enrichplot::pairwise_termsim(ora_go) %>%
  enrichplot::treeplot()

# extract useful gene sets ---------
ora_go@geneSets %>%
  names() %>%
  clusterProfiler::go2term() -> term_collect

read_tsv('../covid19/results/wu-chen-go-pathways-digest.txt') -> digest_pathway

digest_pathway$Term <- reorder(digest_pathway$Term, digest_pathway$order, decreasing = TRUE)

digest_pathway %>%
  rename(count = Up) %>%
  rename(adjust.P = padj.up) %>%
  ggplot(aes(Term, gene_ratio, size = count, color = adjust.P)) +
  geom_point() +
  ggpubr::labs_pubr() +
  coord_flip() +
  expand_limits(y = 0)+
  scale_color_gradient(low = 'red', high = 'black')
