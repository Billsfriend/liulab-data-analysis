# 2022.6.13
# re-analyze LZ B cell scRNA data with Seurat
# from 2020 Holmes JEM paper 
library(Seurat)
library(tidyverse) # contain readr, dplyr, stringr, ggplot2
library(cowplot)
library(org.Hs.eg.db)
library(data.table)
library(DropletUtils)

LZUmi <- fread("data/GSM4148375_LZ2_umi_grch38.txt.gz") # directly read .gz is ok

LZUmi[,1:5]

features <- LZUmi[,1]
features <- as_tibble(features)
features <- separate(features, V1, sep = ';', into = c('ensembl','symbol')) # separate() column in tidyr

LZUmi$symbol <- features$symbol

LZUmi <- unique(LZUmi, by = 'symbol')

LZUmi <- rename(LZUmi, Gene = V1)

LZUmi$Gene <- LZUmi$symbol

# to filter out IgG1/IgM-exclusive expressing cells
LZTb <- transpose(LZUmi, keep.names = 'rn', make.names = 'Gene')

LZTb[1:5,1:5]

LZTb <- LZTb[IGHG1 > 0 | IGHM > 0]

igg1barcode <- LZTb[IGHM == 0, rn]

igMbarcode <- LZTb[IGHG1 == 0, rn]

# Seurat take count matrix without feature names, and take feature names in row.names
Srt <- CreateSeuratObject(LZUmi[,-c('Gene','symbol')], row.names = LZUmi$symbol)

head(Srt@assays$RNA@data)

# quality control of auto-generated metadata
# calculate ratio of mitochondria genes
Srt$mitoRatio <- PercentageFeatureSet(object = Srt, pattern = "^MT-")
Srt$mitoRatio <- Srt@meta.data$mitoRatio / 100

# extract auto-generated data
orimeta <- Srt@meta.data

# calculate and rename some QC data
orimeta$log10GenesPerUMI <- log10(orimeta$nFeature_RNA) / log10(orimeta$nCount_RNA)
orimeta <- orimeta %>%
  dplyr::rename(nUMI = nCount_RNA,
                nGene = nFeature_RNA)
# delete column by dplyr
orimeta <- dplyr::select(orimeta, -orig.ident)

orimeta$cells <- rownames(orimeta)

orimeta$exClass <- NA

orimeta$exClass[which(orimeta$cells %in% igg1barcode)] <- 'IgG1 cells'

orimeta$exClass[which(orimeta$cells %in% igMbarcode)] <- 'IgM cells'

# setting NA as ident is dangerous. set 'NA'.
orimeta$exClass[which(is.na(orimeta$exClass))] <- 'NA'

Srt@meta.data <- orimeta

# write10xCounts('data/AscLZ-B', Srt@assays[["RNA"]]@counts)

Srt <- SetIdent(Srt, value = 'exClass')

write_rds(Srt, 'data/HolmesLZ.rds')


# raw count data cannot plot violin
# VlnPlot(Srt,
#         features = 'IGHG1')
# 
# VlnPlot(Srt,
#         features = 'IGHM')
