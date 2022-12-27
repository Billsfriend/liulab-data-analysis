# 2022.6.13
# re-analyze GC B cell scRNA data with Seurat
# from 2020 Holmes JEM paper 
library(Seurat)
library(tidyverse) # contain readr, dplyr, stringr, ggplot2
library(cowplot)
library(org.Hs.eg.db)
library(data.table)
library(DropletUtils)

gcUmi <- fread("data/GSM4148371_GC2_umi_grch38.txt.gz") # directly read .gz is ok

gcUmi[,1:5]

features <- gcUmi[,1]
features <- as_tibble(features)
features <- separate(features, V1, sep = ';', into = c('ensembl','symbol')) # separate() column in tidyr

gcUmi$symbol <- features$symbol

gcUmi <- unique(gcUmi, by = 'symbol')

gcUmi <- rename(gcUmi, Gene = V1)

gcUmi$Gene <- gcUmi$symbol

gcTb <- transpose(gcUmi, keep.names = 'rn', make.names = 'Gene')

gcTb[1:5,1:5]

gcTb <- gcTb[IGHG1 > 0 | IGHM > 0]

igg1barcode <- gcTb[IGHM == 0, rn]

igMbarcode <- gcTb[IGHG1 == 0, rn]

gcUmi <- gcUmi[,-c('V1','symbol')]

# Seurat take count matrix without feature names, and take feature names in row.names
Srt <- CreateSeuratObject(gcUmi[,-c('V1','symbol')], row.names = gcUmi$symbol)

head(Srt@assays$RNA@data)


# quality control of auto-generated metadata
# calculate ratio of mitochondria genes
Srt$mitoRatio <- PercentageFeatureSet(object = Srt, pattern = "^MT-")
Srt$mitoRatio <- Srt@meta.data$mitoRatio / 100

# extract auto-generated data
orimeta <- Srt@meta.data
orimeta$log10GenesPerUMI <- log10(orimeta$nFeature_RNA) / log10(orimeta$nCount_RNA)
orimeta <- orimeta %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# delete column by dplyr
orimeta <- dplyr::select(orimeta, -orig.ident)

orimeta <- rownames_to_column(orimeta, var = 'cells')

orimeta$exClass <- NA

orimeta$exClass[which(orimeta$cells %in% igg1barcode)] <- 'IgG1 cells'

orimeta$exClass[which(orimeta$cells %in% igMbarcode)] <- 'IgM cells'

# setting NA as ident is dangerous. set 'NA'.
orimeta$exClass[which(is.na(orimeta$exClass))] <- 'NA'

Srt@meta.data <- orimeta

write10xCounts('data/AscGC-B', Srt@assays[["RNA"]]@counts)

write_rds(Srt, 'data/HolmesGC.rds')

VlnPlot(Srt,
        features = 'IGHG1')

VlnPlot(Srt,
        features = 'IGHM')
