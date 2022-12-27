# 2022.3.23
# from 2018 Zhang Cell paper 
library(Seurat)
library(tidyverse)

cnts <- data.table::fread("CRC-I/data/nature/GSE108989_CRC.TCell.S11138.count.txt.gz")

cnts[1:5,1:5]
# find cells belong to our genotyped patients
data <- select(cnts, matches("symbol|0215|0701|0909|1012"))
# SeuratObject input need feature as rows and sample as columns
# Seurat do not allow duplicate feature names
data <- distinct(data, symbol, .keep_all = TRUE)

data %>%
  filter(!is.na(symbol)) %>%
  column_to_rownames('symbol') %>%
  CreateSeuratObject(names.delim = '-',
                     min.cells = 1,
                     min.features = 1) ->
  sobj

sobj

Srt <- CreateSeuratObject(data[,-1], row.names = t(data[, 1]), meta.data = meta)
head(Srt@assays$RNA@data)

# load metadata
meta <- read_delim("CRC-I/data/nature/Tmeta.csv")

# assign I232T genotype in meta
meta %>%
  filter(Patient_ID %in% c("P0215","P0701","P0909","P1012")) %>%
  mutate(genotype = case_when(
  str_detect(Patient_ID, 'P0701|P1012') ~ 'II',
  Patient_ID == 'P0909' ~ 'IT',
  Patient_ID == 'P0215' ~ 'TT'
)) %>%
  mutate(tissue = case_when(
    str_starts(sampleType, 'T') ~ 'Tumor',
    str_starts(sampleType, 'P') ~ 'Blood',
    str_starts(sampleType, 'N') ~ 'Normal'
  ))->
  meta

# set rownames of metadata
meta <- as.data.frame(meta)
rownames(meta) <- meta$UniqueCell_ID

# quality control of auto-generated metadata
# calculate ratio of mitochondria genes
sobj$mitoRatio <- PercentageFeatureSet(object = sobj, pattern = "^MT-")

sobj <- subset(sobj, mitoRatio < 10)

# extract auto-generated data
orimeta <- Srt@meta.data
orimeta$log10GenesPerUMI <- log10(orimeta$nFeature_RNA) / log10(orimeta$nCount_RNA)
orimeta <- orimeta %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# delete column by dplyr
orimeta <- dplyr::select(orimeta, -orig.ident)
# merge auto-generate and attached metadata
mergemeta <- merge(meta, orimeta, by=0)
rownames(mergemeta) <- mergemeta$Row.names
mergemeta <- dplyr::select(mergemeta, -Row.names)
Srt@meta.data <- mergemeta

Srt <- SetIdent(Srt, value = 'ITgeno')
# save genotyped data
write_rds(Srt, 'data/nature/NatureCount.rds')

fwrite(meta[c('majorCluster', 'tissue', 'ITgeno', 'UniqueCell_ID')], 'data/nature/TcellMeta.csv')

# Visualize the number of cell counts per sample
Srt@meta.data %>% 
  ggplot(aes(x=ITgeno, fill=ITgeno)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# ITgeno proprotion in every cell subsets
Srt@meta.data %>% 
  ggplot(aes(x=majorCluster, fill=ITgeno)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

VlnPlot(Srt,
        features = 'FCGR2B',
        group.by = 'Patient_ID',
        split.by = 'ITgeno')

VlnPlot(Srt,
        features = 'FCGR2B', 
        group.by = 'majorCluster',
        split.by = 'ITgeno')
