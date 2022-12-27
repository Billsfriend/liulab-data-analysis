library(tidyverse)
library(SummarizedExperiment)
library(survival)
library(survminer)

coadData <- read_rds('data/COADdata.rds')
rectData <- read_rds('data/READdata.rds')

cluster_signature <- read_delim('ref/CRC_tumor_enrich.txt')

coadFpkm <- as_tibble(coadData@assays@data$fpkm_uq_unstrand)

# there are cases in COAD with multiple samples
colData(coadData) %>%
  as_tibble() %>%
  dplyr::select(patient,
                barcode) %>%
  distinct(patient, .keep_all = TRUE) ->
  coadCases

# keep only one record of multi-case = 456 case
coadFpkm %>%
  select(coadCases$barcode) ->
  coadFpkm

rectFpkm <- as_tibble(rectData@assays@data$fpkm_uq_unstrand)

crcFpkm <- bind_cols(rectFpkm, coadFpkm) 

sample_list <- colnames(crcFpkm)
cell_type_enrichment <- tibble(sample = sample_list)

# log-transformation
crcFpkm <- log(crcFpkm+1)

crcFpkm$symbol <- coadData@rowRanges$gene_name

# save log-fpkm data
write_csv(crcFpkm, 'data/other/tcga-crc-logFpkm.csv.gz')

# calculate enrichment score for each clusters --------
for (i in 1:22) {
  
  clusGenes <- str_subset(cluster_signature[[i]], '[:alnum:]')
  
  cluster_name <- names(cluster_signature[i])
  
  scaleMean <- crcFpkm %>%
    filter(symbol %in% clusGenes) %>% # retain only signature genes
    dplyr::select(-symbol) %>% # remove symbol column
    t() %>% # transpose to make gene as column
    scale() %>% # scale at column (gene) level 
    t() %>% # transpose to make sample column
    as_tibble() %>%
    summarize_all(mean) %>% # compute mean for each sample
    as.numeric() %>%
    as_tibble()
  
  colnames(scaleMean) <- cluster_name
  
  cell_type_enrichment %>%
    bind_cols(scaleMean) ->
    cell_type_enrichment
}

cell_type_enrichment %>%
  mutate(patient = str_trunc(sample, 12, ellipsis = '')) ->
  cell_type_enrichment

cell_type_enrichment[-9] ->
  cell_type_enrichment

write_csv(cell_type_enrichment, 'results/tcga_coad+read_cell_type_enrichment.csv')

cell_type_enrichment <- read_csv('results/tcga_coad+read_cell_type_enrichment.csv')

# use specific cutoff to subset low-high enrichment -----------
# lower 25% as low-enrichment, upper 75% as high-enrichment
# or < mean - 0.75SD as low, > mean + 0.75SD as high
quandral_enrich <- function(sample_vector,
                            mean_vector,
                            name,
                            mode = 'sd') {
  data <- tibble(sample = sample_vector,
                 mean = mean_vector)
  
  if (mode == 'quantile') {
    data %>%
      summarise(x = quantile(mean, probs = c(0.25, 0.75))) ->
      quant_cluster
  } else {
    quant_cluster = tibble(x = c(0.75, -0.75))
  }
  
  data %>%
    filter(mean >= max(quant_cluster)) %>%
    mutate(enrich = 'high') ->
    highEnrich
  
  data %>%
    filter(mean <= min(quant_cluster)) %>%
    mutate(enrich = 'low') ->
    lowEnrich
  
  data %>%
    filter(mean > min(quant_cluster) & mean < max(quant_cluster)) %>%
    bind_rows(lowEnrich) %>%
    bind_rows(highEnrich) ->
    scale_cluster
  
  names(scale_cluster) <- c('sample', 'mean', name)
  
  scale_cluster <- left_join(data, scale_cluster)
  
  scale_cluster[-2]
}

clusterEnrich <- tibble(sample = cell_type_enrichment$patient)

for (i in 1:21) {
  t <- quandral_enrich(
    cell_type_enrichment$patient,
    cell_type_enrichment[[i]],
    str_sub(names(cell_type_enrichment[i]),
            start = 6L),
    mode = 'quantile'
  )
  
  clusterEnrich <- left_join(clusterEnrich, t)
}

clusterEnrich <- left_join(clusterEnrich, cell_type_enrichment, by = c('sample' = 'patient'))

clusterEnrich %>%
  colnames() %>%
  str_replace('-','_') ->
  colnames(clusterEnrich)

write.csv(clusterEnrich, 'results/enrich_25percent_cutoff.csv')

# visualize the enrichment result -------
clusterEnrich %>%
  ggplot(aes(x = CD4_TCF7, y = hT04_CD4_TCF7)) +
  geom_boxplot(color = hT04_CD4_TCF7)

# merge information of two cluster
scaleSpp1 <- scaleSpp1[-5]

scaleClus <- left_join(scaleSpp1, scaleC1qc)

scaleClus <- scaleClus %>%
  dplyr::select(patient, C1QCcluster, SPP1cluster)

# extract survival data ----------
colData(coadData) %>%
  as_tibble() %>%
  dplyr::select(patient,
         days_to_last_follow_up,
         days_to_death,
         vital_status) %>%
  distinct(patient, .keep_all = TRUE) ->
  coadSurvival

colData(rectData) %>%
  as_tibble() %>%
  dplyr::select(patient,
                days_to_last_follow_up,
                days_to_death,
                vital_status) %>%
  distinct(patient, .keep_all = TRUE) ->
  rectSurvival

crcSurvival <- bind_rows(coadSurvival, rectSurvival)

# get a column of follow-up days
crcSurvival %>%
  filter(!is.na(days_to_death)) %>%
  mutate(days_to_last_follow_up = days_to_death) %>%
  bind_rows(filter(crcSurvival, is.na(days_to_death))) ->
  crcSurvival

crcSurvival <- crcSurvival %>%
  mutate(dead = vital_status == 'Dead') %>%
  select(-days_to_death)

write_csv(crcSurvival,'results/coad+read_survival.csv')

crcSurvival <- read_csv('results/coad+read_survival.csv')

crcSurvCluster <- crcSurvival %>%
  left_join(clusterEnrich, by = c('patient' = 'sample'))

crcSurvival[1:4] -> crcSurvival
