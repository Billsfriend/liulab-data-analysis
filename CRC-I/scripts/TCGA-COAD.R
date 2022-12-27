library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)

# browseVignettes('TCGAbiolinks')

# query and download data ----------
# Colon Adenocarcinoma from TCGA
query <- GDCquery(project = 'TCGA-COAD',
                  data.category = 'Transcriptome Profiling',
                  data.type = 'Gene Expression Quantification',
                  sample.type = 'Primary Tumor')

# api download can be interrupted frequently
repeat try(GDCdownload(query, files.per.chunk = 5))

COADdata <- GDCprepare(query)

# remove single-stranded and quasi-normalized FPKM assays
assays(COADdata) <- assays(COADdata)[6]

clinical <- COADdata@colData %>%
  as.data.frame() %>%
  as_tibble()

survival <- clinical %>%
  select(c(patient, days_to_last_follow_up, vital_status)) %>%
  filter(!is.na(days_to_last_follow_up)) %>%
  filter(!is.na(vital_status))

write_rds(COADdata, 'data/COADdata.rds')

# Rectum Adenocarcinoma
query <- GDCquery(project = 'TCGA-READ',
                  data.category = 'Transcriptome Profiling',
                  data.type = 'Gene Expression Quantification',
                  sample.type = 'Primary Tumor')

# use official GDC-transfer tool instead. it's stable

READdata <- GDCprepare(query)

# keep only normalized fpkm data
assays(READdata) <- assays(READdata)[6]

write_rds(READdata, 'data/READdata.rds')

# CPTAC(clinical proteomic tumor analysis consortium)-3 scRNA-seq data
query3 <- GDCquery(project = 'CPTAC-3',
                  data.category = 'Transcriptome Profiling',
                  data.format = 'hdf5')

# will throw error
# GDCdownload(query3, files.per.chunk = 1)

# not preparing loom files
# cptac3 <- GDCprepare(query3)

# install.packages('hdf5r')
# remotes::install_github("mojaveazure/seurat-disk")
library(Seurat)
library(SeuratDisk)

loom <- Connect(filename = '~/Documents/DataAnalysis/CRC-I/GDCdata/CPTAC-3/harmonized/Transcriptome_Profiling/Single_Cell_Analysis/052a06ca-73ea-4a40-aa6a-7118deee3078/seurat.loom', mode = 'r')

loomSeurat <- as.Seurat(loom)

ensemble <- loomSeurat@assays[["RNA"]]@counts@Dimnames[["features"]]

library(biomaRt)


VlnPlot(loomSeurat, 'SPP1')
