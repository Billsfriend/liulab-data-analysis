# for CRC0926 human sample
library(Seurat)
library(tidyverse)
library(tidyseurat)

data.table::fread('CRC-I/data/seekgene/CRC_TT_0926.gz') -> data

data[1:5,1:5]

data %>% 
  column_to_rownames('V1') %>%
  t() %>%
  CreateSeuratObject(min.cells = 1,
                     min.features = 1) ->
  sobj

sobj$mitoRatio <- PercentageFeatureSet(object = sobj, pattern = "^MT-")

sobj %>%
  filter(nFeature_RNA >= 200 & mitoRatio < 10) ->
  sobj

write_rds(sobj, 'CRC-I/data/seekgene/crc0926.rds','gz')

sobj@meta.data -> meta

meta %>%
  rownames_to_column('barcode') %>%
  filter(barcode == 'AATGCCAATGGAGGCA')
