# remotes::install_github("MarioniLab/DropletUtils")
library(DropletUtils)
library(Seurat)
library(data.table)

NatureCount <- readRDS("~/Documents/DataAnalysis/CRC-I/data/nature/NatureCount.rds")

# need a non-exist directory
write10xCounts('~/Documents/DataAnalysis/CRC-I/data/nature/10xCount', NatureCount@assays[["RNA"]]@counts)
