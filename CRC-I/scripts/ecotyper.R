library(Seurat)
library(tidyverse)
library(ggpubr)

sgfiles <- list.files('ecotyper-master/RecoveryOutput/crc0926_immune_mat', recursive = TRUE, full.names = TRUE, pattern = 'state_assignment.txt') %>%
  map(read_tsv) %>%
  purrr::reduce(bind_rows)

zyfiles <- list.files('ecotyper-master/RecoveryOutput/zy2020_tme_mat', recursive = TRUE, full.names = TRUE, pattern = 'state_assignment.txt') %>%
  map(read_tsv) %>%
  purrr::reduce(bind_rows)

sgfiles <- read_tsv('CRC-I/results/crc0926_immune_eco_anno.tsv') %>%
  right_join(sgfiles)

eco_data <- read_tsv('CRC-I/results/zy2020_immune_eco_anno.tsv') %>%
  right_join(zyfiles) %>%
  bind_rows(sgfiles)

eco_data <- eco_data %>%
  mutate(genotype = case_when(
    genotype == 'TT' ~ 'IT',
    TRUE ~ genotype
  )) 

eco_data %>%
  filter(CellType != 'PMNs' & CellType != 'Mast.cells') %>%
  ggplot(aes(genotype, fill = CellType)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Paired') +
  theme_pubr() +
  ylab('Proportion in TME immune cells')

eco_data %>%
  filter(CellType != 'PMNs' & CellType != 'Mast.cells') %>%
  ggplot(aes(genotype, fill = State)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  facet_wrap(~ CellType) +
  scale_fill_brewer(type = 'qual', palette = 'Paired') +
  theme_pubr() +
  ylab('Proportion in TME immune cells')

eco_data %>%
  filter(str_detect(CellType, 'B.cell')) %>%
  ggplot(aes(genotype, fill = State)) +
  geom_bar(position = 'fill') +
  coord_flip() +
  scale_fill_brewer(type = 'qual', palette = 'Paired') +
  theme_pubr() +
  labs_pubr()+
  ylab('Proportion in all TME plasma cells')
