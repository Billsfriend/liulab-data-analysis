# visualize_cell_fraction.R

library(Seurat)
library(tidyverse)
library(ggpubr)
library(pheatmap)

read_rds('CRC-I/data/crc0926_singler.rds') -> sg_sobj
read_rds('CRC-I/data/Zhang-Yu-2020/tenx_singler.rds') -> zy_sobj
read_rds('CRC-I/data/scMerge_CRC221230.rds') -> scM_sobj

sg_sobj@meta.data %>%
  rownames_to_column('rownames') %>%
  bind_rows(rownames_to_column(zy_sobj@meta.data, 'rownames')) %>%
  mutate(genotype = case_when(batch == 'seekgene' ~ 'TT', TRUE ~ genotype)) ->
  separate_meta

separate_meta %>%
  ggplot(aes(x = genotype, fill = singler_labels))+
  geom_bar()+
  coord_flip()+
  ggpubr()+
  ggtitle('TME immune cells in CRC from different FCGR2B-I232T genotypes')

separate_meta %>%
  ggplot(aes(x = genotype, fill = Sub_Cluster))+
  geom_bar()+
  coord_flip()+
  ggpubr()
