library(Seurat)
library(ggpubr)
library(tidyverse)
library(readxl)

# read from excel ------------
meta_from_paper <- excel_sheets('data/source_ESM.xlsx') %>%
  set_names() %>%
  map(read_excel,
      path = 'data/source_ESM.xlsx',
      skip = 1)

pub_cell_meta <- meta_from_paper[["Metadata for fibroblasts"]] %>%
  select(c('cell', 'celltype', 'tissue', 'sample')) %>%
  mutate(cell = str_replace_all(cell, '\\.', '-'))

sobj <- read_rds('data/zhang_B_cell.rds')

# remove 'T' from 'P130T...' etc
meta <- sobj@meta.data %>%
  mutate(barcode = str_replace(barcode, 'T', ''))

meta_t_new <- left_join(meta, pub_cell_meta, by = c('barcode' = 'cell')) 

rownames(meta_t_new) <- rownames(meta)

celltype_na <- filter(meta_t_new, is.na(celltype)) %>%
  select(-c(celltype, tSNE1, tSNE2))

celltype_na$barcode <- rownames(celltype_na)

celltype_na_new <- left_join(celltype_na,
                             pub_cell_meta,
                             by = c('barcode' = 'cell'))

meta_t_new %>%
  filter(!is.na(celltype)) %>%
  bind_rows(celltype_na_new) ->
  meta_t_newest

rownames(meta_t_newest) <- rownames(meta)

sobj@meta.data <- meta_t_newest

ggplot(meta_t_newest) +
  geom_point(aes(x = tSNE1, y = tSNE2, color = celltype, group = genotype), size = 0.5) +
  theme_classic()

ggplot(meta_t_newest) +
  geom_bar(aes(x = genotype, fill = celltype)) +
  scale_fill_brewer(palette = 'Paired') -> p1

ggplot(meta_t_newest) +
  geom_bar(aes(x = genotype, fill = celltype), position = 'fill') +
  scale_fill_brewer(palette = 'Paired') +
  labs(y = 'Proportion') -> p2

p1+p2
write_csv(meta_t_newest, 'results/zhang_B_cell_meta.csv')

write_rds(sobj, 'data/zhang_myeloid.rds')

# dry run on CD45- cells
ITlist <- read.table('data/ITlist.txt')

ITlist %>%
  mutate(V1 = case_when(
    V1 == 'P128T' ~ V1,
    TRUE ~ str_remove(V1, 'T')
  )) -> ITlist

metadata <- separate(pub_cell_meta,
                     col = cell,
                     into = c('patient', 'letter', 'barcode'), 
                     remove = FALSE)

table(metadata$patient)

metadata %>%
  filter(!str_detect(patient, 'N')) %>%
  mutate(genotype = case_when(
    patient %in% c('P8', 'P63', 'P130T') ~ 'TT',
    patient %in% ITlist$V1 ~ 'IT',
    TRUE ~ 'II'
  )) -> metadata

ggplot(metadata) -> p
p + geom_bar(aes(x = genotype, fill = celltype, group = celltype)) -> p1

p + geom_bar(aes(x = genotype, fill = celltype, group = celltype), position = 'fill') +
  labs(y = 'proportion') -> p2

p1 + scale_fill_brewer(palette = 'Paired') + p2 + scale_fill_brewer(palette = 'Paired')

# statistics for all CD45+ cells --------
cd45_cell_meta <- meta_from_paper[["Metadata for fibroblasts"]] %>%
  select(c('cell', 'celltype', 'tissue', 'sample')) %>%
  filter(tissue == 'Tumor')

ITlist %>%
  mutate(noT = str_replace(V1, 'T', '')) ->
  ITlist

cd45_cell_meta$genotype <- 'II'

cd45_cell_meta$genotype[which(cd45_cell_meta$sample %in% c('P8', 'P63', 'P130'))] <- 'TT'

cd45_cell_meta$genotype[which(cd45_cell_meta$sample %in% ITlist$V2)] <- 'IT'

cd45_cell_meta %>%
  ggplot() +
  geom_bar(aes(x = genotype, fill = celltype, group = celltype)) -> p1

cd45_cell_meta %>%
  ggplot() +
  geom_bar(aes(x = genotype, fill = celltype, group = celltype), position = 'fill') +
  labs(y = 'proportion') -> p2

# patchwork is so easy!
p1 + p2

write_csv(cd45_cell_meta, 'results/zhang_CD45+_cell_meta.csv')

# epithelial program vs genotype -----
meta_from_paper$`heterogeneity of epithelial` %>%
  mutate(genotype = case_when(
    sample %in% c('P8', 'P63', 'P130') ~ 'TT',
    sample %in% ITlist$noT ~ 'IT',
    TRUE ~ 'II'
  )) -> epi_hetero_score

ggplot(epi_hetero_score, aes(genotype, `Inter-tumor_heterogeneity`)) +
  geom_jitter(aes(color = genotype),width = 0.01) +
  stat_summary(geom = 'crossbar', fun.data = 'mean_cl_boot', width = 0.1) +
  ggtitle('Epithelial heterogeneity relative to average ESCC') +
  stat_compare_means(comparisons = list(
    c('IT','TT'),c('II','TT')))

meta_from_paper$`Fig. 3c` -> prog_pathways
meta_from_paper$`Fig. 3e; Fig. 3f` -> prog_by_sample
  

prog_by_sample %>%
  mutate(genotype = case_when(
    Sample %in% c('P8', 'P63', 'P130') ~ 'TT',
    Sample %in% ITlist$noT ~ 'IT',
    TRUE ~ 'II'
  )) -> prog_by_sample

ggplot(prog_by_sample, aes(x = genotype, y =`Program score`)) +
  geom_jitter(aes(color = genotype),width = 0.1) +
  stat_summary(geom = 'crossbar', color = 'red', fun = 'mean', width = 0.2) +
  facet_wrap(vars(Program), scales = 'free')

# Treg score by samples --------
meta_from_paper$`Fig. 4b` %>%
  mutate(genotype = case_when(
    Sample %in% c('P8', 'P63', 'P130') ~ 'TT',
    Sample %in% ITlist$noT ~ 'IT',
    TRUE ~ 'II'
  )) %>%
  ggplot(aes(genotype, Treg_score)) +
  geom_jitter(aes(color = genotype),width = 0.1) +
  stat_summary(geom = 'crossbar', color = 'red', fun = 'mean', width = 0.2)+
  stat_compare_means(comparisons = list(
    c('IT','TT'),c('II','TT'))) +
  ggtitle('Treg score of TME T cells')

# PD-L1/2 expression in tDC
meta_from_paper$`Fig. 4h` %>%
  filter(Tissue == 'Tumor') %>%
  mutate(genotype = case_when(
    Sample %in% c('P8', 'P63', 'P130') ~ 'TT',
    Sample %in% ITlist$noT ~ 'IT',
    TRUE ~ 'II'
  )) %>%
  ggplot(aes(genotype, Expression_in_tDC)) +
  geom_jitter(aes(color = genotype),width = 0.1) +
  stat_summary(geom = 'crossbar', color = 'red', fun = 'mean', width = 0.2)+
  stat_compare_means(comparisons = list(
    c('IT','TT'),c('II','TT'))) +
  ggtitle('PD-L1/2 expression of TME tDC cells')+
  facet_wrap(vars(Gene))


meta_from_paper$`Supplementary Fig. 2h` -> sup2h

first(sup2h) %>% as.character() -> colnames(sup2h)

slice(sup2h, -1) -> sup2h

sup2h %>%
  separate(cell, into = c('sample','letter','barcode'), remove = FALSE) %>%
  mutate(genotype = case_when(
    sample %in% c('P8', 'P63', 'P130T') ~ 'TT',
    sample %in% ITlist$V1 ~ 'IT',
    TRUE ~ 'II'
  )) -> sup2h

ggplot(sup2h, aes(genotype, fill = Epi1)) +
  scale_fill_manual(name = 'epithelial cells', label = c('Not Epi1 program', 'Epi1 program'), values = c('grey', 'red')) -> g

g + geom_bar(position = 'fill') + ylab(label = 'proportion') -> g1
g + geom_bar() -> g2
g2 + g1
