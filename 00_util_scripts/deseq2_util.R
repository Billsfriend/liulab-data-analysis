library(DESeq2)
library(pheatmap)
library(clusterProfiler)
library(tidyverse)

read_tabulars <- function(i){
  read_tsv(paste0('../covid19/data/zhang-jiang2022inactivated/2dose',i,'.tabular'))}

read_tabulars(1) %>%
  left_join(read_tabulars(2)) %>%
  left_join(read_tabulars(3)) %>%
  left_join(read_tabulars(4)) %>%
  left_join(read_tabulars(5)) %>%
  left_join(read_tabulars(6)) -> count_data

colnames(count_data) <- c('geneid', 'sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6')

bitr(count_data$geneid,
     fromType = 'ENTREZID',
     toType = 'SYMBOL',
     OrgDb = 'org.Hs.eg.db') -> gene_names

genotype <- factor(c('GG','RR','GG','GG','GG','GG'))

FCGR2B <- factor(c('II','II','II','IT','II','II'))

col_data <- data.frame(genotype, FCGR2B, row.names = c('sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'))

count_data %>%
  column_to_rownames('geneid') %>%
  DESeqDataSetFromMatrix(colData = col_data,
                       design = ~ FCGR2B + genotype) -> dds

# pre-filtering data of too low counts
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

# relevel to specify the control group (reference)
levels(dds$genotype) <- c('RR', 'GG')
dds$genotype <- relevel(dds$genotype, ref = "GG")

rld <- rlog(object = dds, blind=FALSE)
vsd <- vst(object=dds, blind=TRUE)

dds <- DESeq(dds)
res <- results(dds)
summary(res)

normalized_counts <- assay(rld)

as.data.frame(res) %>%
  rownames_to_column('ENTREZID') %>%
  left_join(gene_names) -> deg_tibble

write_csv(deg_tibble, '../covid19/results/zhang-jiang-inactivated-deseq2.csv')

deg_tibble %>%
  mutate(type = case_when(
    log2FoldChange > 1.5 & padj < 0.1 ~ 'Enriched in RR',
    log2FoldChange < -1.5 & padj < 0.1 ~ 'Depleted in RR',
    TRUE ~ 'NS'
  )) %>%
  mutate(sig_symbol = case_when(
    type != 'NS' ~ SYMBOL,
    TRUE ~ ''
  )) -> deg_tibble

ggplot(deg_tibble, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(color = type)) +
  geom_vline(xintercept = c(1.5, -1.5), linetype = 'dashed') +
  geom_hline(yintercept = 3.5, linetype = 'dashed') +
  scale_color_manual(values = c('red', 'blue', 'grey'))

# clusterProfiler -----------
deg_tibble %>%
  filter(log2FoldChange > 0, padj < 0.1) %>%
  pull(SYMBOL) -> up_list

deg_tibble %>%
  filter(log2FoldChange < 0, padj < 0.1) %>%
  pull(SYMBOL) -> down_list

enrichGO(down_list,
         OrgDb = 'org.Hs.eg.db',
         keyType = 'SYMBOL',
         ont = 'ALL',
         universe = deg_tibble$SYMBOL,
         readable = TRUE) -> ora_go

dotplot(ora_go)

cnetplot(ora_go)

enrichplot::pairwise_termsim(ora_go) %>%
  enrichplot::treeplot()
