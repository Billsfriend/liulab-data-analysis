# based on Seurat official vignette
# https://satijalab.org/seurat/articles/de_vignette.html#perform-default-differential-expression-tests-1

install.packages("devtools")
devtools::install_github('satijalab/seurat-data')

library(Seurat)
library(SeuratData)

InstallData('pbmc3k')
pbmc <- LoadData("pbmc3k", type = "pbmc3k.final")

# list options for groups to perform differential expression on
levels(pbmc)

# Find differentially expressed features between CD14+ and FCGR3A+ Monocytes
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono")
# view results
head(monocyte.de.markers)

# Find differentially expressed features between CD14+ Monocytes and all other cells, only
# search for positive markers
monocyte.de.markers <- FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = NULL, only.pos = TRUE)
# view results
head(monocyte.de.markers)

# Pre-filter features that are detected at <50% frequency in either CD14+ Monocytes or FCGR3A+
# Monocytes
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.pct = 0.5))

# Pre-filter features that have less than a two-fold change between the average expression of
# CD14+ Monocytes vs FCGR3A+ Monocytes
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", logfc.threshold = log(2)))

# Pre-filter features whose detection percentages across the two groups are similar (within
# 0.25)
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", min.diff.pct = 0.25))

# Increasing min.pct, logfc.threshold, and min.diff.pct, will increase the speed of DE
# testing, but could also miss features that are prefiltered

# Subsample each group to a maximum of 200 cells. Can be very useful for large clusters, or
# computationally-intensive DE tests
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", max.cells.per.ident = 200))

# Test for DE features using the DESeq2 package. Throws an error if DESeq2 has not already
# been installed Note that the DESeq2 workflows can be computationally intensive for large
# datasets, but are incompatible with some feature pre-filtering options We therefore suggest
# initially limiting the number of cells used for testing
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "DESeq2", max.cells.per.ident = 50))

# Test for DE features using the MAST package
head(FindMarkers(pbmc, ident.1 = "CD14+ Mono", ident.2 = "FCGR3A+ Mono", test.use = "MAST"))