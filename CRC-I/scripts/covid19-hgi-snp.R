library(data.table)
library(tidyverse)

my_snp <- c('rs1050501')

read.csv('CRC-I/covid10-hgi-my-snp.csv') -> t

# main ----------
# very severe vs population
fread(input = 'https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats/COVID19_HGI_A2_ALL_leave_23andme_20220403.tsv.gz') %>%
  filter(rsid %in% my_snp) -> all_hgi_A2

# hospitalized cov vs non-hospitalized cov
fread(input = 'https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats/COVID19_HGI_B1_ALL_leave_23andme_20220403.tsv.gz') %>%
  filter(rsid %in% my_snp) -> all_hgi_B1

# hospitalized cov vs population
fread(input = 'https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats/COVID19_HGI_B2_ALL_leave_23andme_20220403.tsv.gz') %>%
  filter(rsid %in% my_snp) -> all_hgi_B2

# cov vs population
fread(input = 'https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/main/sumstats/COVID19_HGI_C2_ALL_leave_23andme_20220403.tsv.gz') %>%
  filter(rsid %in% my_snp) -> all_hgi_C2

# population-specific ------------
get_popu_spec_hgi <- function(popu, group = 'B2'){
  fread(input = str_glue('https://storage.googleapis.com/covid19-hg-public/freeze_7/results/20220403/pop_spec/sumstats/COVID19_HGI_{group}_ALL_{popu}_leave23andme_20220403.tsv.gz')) %>%
    filter(rsid %in% my_snp)
}

popu_list <- c('afr', 'amr', 'eas', 'eur', 'sas')

group_list <- c('A2', 'B2', 'C2')

get_popu_spec_hgi('afr', 'A2') -> afr_A2
get_popu_spec_hgi('amr', 'A2') -> amr_A2
get_popu_spec_hgi('eas', 'A2') -> eas_A2
get_popu_spec_hgi('eur', 'A2') -> eur_A2
get_popu_spec_hgi('sas', 'A2') -> sas_A2

expand_grid(popu = popu_list, group = group_list) -> param_list

map2(param_list$popu, param_list$group, get_popu_spec_hgi) -> b2

param_list$result <- b2

param_list %>%
  rowwise() %>%
  filter(size_sum(result) == '[1 × 15]') ->
  param_list

param_list %>%
  pull(result) %>%
  bind_rows() %>%
  bind_cols(param_list[1:2]) ->
  subpopu_hgi

select(all_hgi_C2, !contains('lmso')) %>%
  add_column(popu = 'all', group = 'C2')->
  all_hgi_C2

bind_rows(subpopu_hgi, all_hgi_C2) -> hgi_my_snp

write_csv(hgi_my_snp, 'covid19/results/hgi_FCGR2B_I232T.csv')

hgi_my_snp <- read_csv('covid19/results/hgi_FCGR2B_I232T.csv')

# simple t test ---------
test_popu_proportion <- function(size1, event1, size2, event2){
  array1 <- c(rep.int(1, event1), rep.int(0, size1 - event1))
  array2 <- c(rep.int(1, event2), rep.int(0, size2 - event2))
  test_result <- t.test(array1, array2)
  test_result$p.value
}

read_tsv('1000genome_igg1.txt') -> region_freq

region_freq %>%
  rowwise() %>%
  mutate(vs_CHB_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[1], 
    region_freq$minorCount[1])) %>%
  mutate(vs_CHS_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[2], 
    region_freq$minorCount[2])) %>%
  mutate(vs_CDX_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[3], 
    region_freq$minorCount[3]))%>%
  mutate(vs_BEB_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[4], 
    region_freq$minorCount[4]))%>%
  mutate(vs_JPT_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[5], 
    region_freq$minorCount[5]))%>%
  mutate(vs_KHV_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[6], 
    region_freq$minorCount[6]))%>%
  mutate(vs_PJL_p_val = test_popu_proportion(
    nSize,
    minorCount,
    region_freq$nSize[7], 
    region_freq$minorCount[7])) -> p_val_mat

write_csv(p_val_mat, '1000genome_region_pval.csv')

region_freq$nSize -> size_list
region_freq$minorCount -> minor_count_list

matrix(nrow = 7, ncol = 7) -> p_val_mat

for (i in 1:7) {
  ref_size <- size_list[i]
  ref_minor <- minor_count_list[i]
  for (j in 1:7) {
    test_popu_proportion(ref_size,
                         ref_minor,
                         size_list[j],
                         minor_count_list[j]) ->
      p_val_mat[i, j]
  }
}
