library(SingleCellExperiment)
library(scater)
library(scMerge)
library(tidyverse)

## Subsetted mouse ESC data
data("example_sce", package = "scMerge")

example_sce = runPCA(example_sce, exprs_values = "logcounts")

plotPCA(
  example_sce, 
  colour_by = "cellTypes", 
  shape_by = "batch")

# 1. obtaining negative ctrl: SEG ---------
## single-cell stably expressed gene list
data("segList_ensemblGeneID", package = "scMerge")
head(segList_ensemblGeneID$mouse$mouse_scSEG)

# 2. unsupervised scMerge ----------
# kmeansK value means there are 3 types of cell in each sample
scMerge_unsupervised <- scMerge(
  sce_combine = example_sce, 
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = "scMerge_unsupervised"
  )

scMerge_unsupervised = runPCA(scMerge_unsupervised, exprs_values = "scMerge_unsupervised")
plotPCA(
  scMerge_unsupervised, 
  colour_by = "cellTypes", 
  shape_by = "batch")
# out of RAM computation ----
library(HDF5Array)
library(DelayedArray)

DelayedArray:::set_verbose_block_processing(TRUE) ## To monitor block processing 

hdf5_input = example_sce

assay(hdf5_input, "counts") = as(counts(hdf5_input), "HDF5Array")
assay(hdf5_input, "logcounts") = as(logcounts(hdf5_input), "HDF5Array")
system.time(scMerge_hdf5 <- scMerge(
  sce_combine = hdf5_input, 
  ctl = segList_ensemblGeneID$mouse$mouse_scSEG,
  kmeansK = c(3, 3),
  assay_name = "scMerge_hdf5",
  BSPARAM = BiocSingular::RandomParam(),
  BPPARAM = BiocParallel::MulticoreParam(workers = 6),
  verbose = TRUE)) 
