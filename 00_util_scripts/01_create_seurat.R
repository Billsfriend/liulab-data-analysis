library(Seurat)
library(harmony)
library(singleCellTK)
library(data.table)
library(tidyverse)
library(tidyseurat)

# read in data ---------
# use url when possible
fread('data/nature/GSE108989_CRC.TCell.S11138.count.txt.gz') -> mtx

mtx[1:5,1:5]

mtx %>%
  column_to_rownames('geneID') %>%
  select(-symbol) %>%
  CreateSeuratObject() -> sobj

VlnPlot(sobj, 'nCount_RNA')

mean(sobj$nCount_RNA)

write_rds()

read_rds('../covid19/data/Final_nCoV_0716_upload_BGI.RDS') -> sobj

glimpse(sobj@meta.data)

table(sobj$batch)

str_which(sobj$batch, 'COV') -> cov_cells

subset(sobj, cell = cov_cells) -> sobj

mutate(sobj, genotype = case_when(
  str_detect(batch, 'COV-4') ~ 'RR',
  TRUE ~ 'GG'
)) %>%
  mutate(severity = case_when(
    str_detect(batch, 'COV-5') ~ 'severe',
               TRUE ~ 'mild')
  ) -> sobj

sobj %>%
  separate(col = 'batch', into = c('group', 'patient', 'time'), remove = FALSE) -> sobj

write_rds(sobj, '../covid19/data/bgi-cov-geno.rds')

sobj <- read_rds('../covid19/data/bgi-cov-geno.rds')

subset(sobj, severity == 'mild') -> sobj

ggplot(sobj, aes(Stage, batch)) +
  geom_count()

filter(sobj, time == 'D16') %>%
  ggplot(aes(genotype, fill = cell_type)) +
  scale_fill_viridis_d(option = 'turbo') -> p

p + geom_bar() + NoLegend() -> p1
p + geom_bar(position = 'fill') + ylab('proportion') -> p2

p1 + p2

# plot every cell types in facets
sobj %>%
  filter(time == 'D16') %>%
  group_by(genotype, cell_type) %>%
  tally() -> count_data

count_data %>%
  group_by(genotype) %>%
  summarise(sum(n)) -> sum_genotype

count_data %>%
  rowwise() %>%
  mutate(percent = case_when(
    genotype == 'GG' ~ n/sum_genotype$`sum(n)`[1],
    TRUE ~ n/sum_genotype$`sum(n)`[2]
  )) -> percent_data

ggplot(percent_data, aes(genotype, percent))+
  geom_col(aes(fill = genotype), position = 'dodge') +
  #coord_flip() +
  facet_wrap(vars(cell_type), scales = 'free')

# gene-level analysis ---------
SetIdent(sobj, value = 'genotype') -> sobj

SplitObject(sobj, 'cell_type') -> cell_type_list

enrichR::listEnrichrDbs() -> t

get_path_enrich <- function(fobj){
  Seurat::DEenrichRPlot(fobj,
                        ident.1 = 'RR',
                        ident.2 = 'GG',
                        enrich.database = 'GO_Biological_Process_2021',
                        max.genes = 10000,
                        num.pathway = 24,
                        return.gene.list = TRUE)
}

map(cell_type_list, safely(get_path_enrich)) %>%
  transpose() %>%
  pluck('result') %>%
  compact()-> enrich_result

names(enrich_result) -> enriched_types

transpose(enrich_result) %>%
  pluck('pos') %>%
  set_names(enriched_types) %>%
  compact() -> pos_enrich_result

transpose(enrich_result) %>%
  pluck('neg') %>%
  set_names(enriched_types) %>%
  compact() -> neg_enrich_result

write_rds(enrich_result, '../covid19/results/zhu-liu-enrich-result.rds')

read_rds('../covid19/results/zhu-liu-enrich-result.rds') -> enrich_result
