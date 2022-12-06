# https://rpubs.com/tiagochst/TCGAworkshop is obsolete
# BiocManager::install('TCGAbiolinks')
# BiocManager::install('MultiAssayExperiment')
# BiocManager::install('maftools')
# BiocManager::install('ComplexHeatmap')
library(TCGAbiolinks)
# browseVignettes('TCGAbiolinks')
library(MultiAssayExperiment)
library(maftools)
library(dplyr)
library(ComplexHeatmap)
library(data.table)
library(DT) # interactive table

# Access indexed clinical data
clinical <- GDCquery_clinic("TCGA-COAD")
head(clinical)

# extract one patient
clinical %>% 
  dplyr::filter(submitter_id == "TCGA-AA-3562") %>% 
  t %>% 
  as.data.frame

# or, same data can be download as tsv from BCR biotab
query <- GDCquery(project = "TCGA-ACC", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)

# download two sample of RNA-seq counts
query.exp.hg38 <- GDCquery(project = "TCGA-GBM", 
                           data.category = "Transcriptome Profiling", 
                           data.type = "Gene Expression Quantification", 
                           workflow.type = "STAR - Counts",
                           barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
GDCdownload(query.exp.hg38)
raw.counts <- GDCprepare(query = query.exp.hg38, summarizedExperiment = FALSE)
head(raw.counts)

# or download normalized FPKM data
# ...NO LONGER PROVIDED

# access the harmonized database (legacy = FALSE) and search for all DNA methylation data for recurrent glioblastoma multiform (GBM) and low grade gliomas (LGG) samples
query <- GDCquery(
  project = c("TCGA-GBM", "TCGA-LGG"),
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  sample.type = "Recurrent Tumor"
)
data.table(
  getResults(query), 
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)

# DNA-methylation & RNA-seq data from TCGA-COAD
query.met <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450")
)
query.exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts"
)

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(
  substr(getResults(query.met, cols = "cases"), 1, 12),
  substr(getResults(query.exp, cols = "cases"), 1, 12)
)

# Only seelct the first 5 patients
query.met <- GDCquery(
  project = "TCGA-COAD",
  data.category = "DNA Methylation",
  legacy = FALSE,
  platform = c("Illumina Human Methylation 450"),
  barcode = common.patients[1:5]
)

query.exp <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification", 
  workflow.type = "STAR - Counts",
  barcode = common.patients[1:5]
)

# retrieve data by getResults()
datatable(
  getResults(query.met, cols = c("data_type","cases")),
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
  rownames = FALSE
)
datatable(
  getResults(query.exp)
)

# download from api/client by GDCdownload
GDCdownload(query.exp)

# after downloading, GDCprepare can produce summarizedExperiment by default
COADexp <- GDCprepare(query.exp)

# proteome data is available too, in some project
query.rppa <- GDCquery(
  project = "TCGA-ESCA", 
  data.category = "Proteome Profiling",
  data.type = "Protein Expression Quantification"
)
GDCdownload(query.rppa) 
rppa <- GDCprepare(query.rppa)

# scRNA-seq data is only available in CPTAC-3 project
query.raw.counts <- GDCquery(
  project = "CPTAC-3", 
  data.category = "Transcriptome Profiling",
  legacy = FALSE,
  access = "open",
  data.type = "Gene Expression Quantification",
  barcode = c("CPT0167860015","CPT0206880004"),
  workflow.type = "CellRanger - 10x Raw Counts"
)  
GDCdownload(query.raw.counts)
raw.counts.list <- GDCprepare(query.raw.counts)


# download SNV data
query <- GDCquery(
  project = "TCGA-COAD", 
  data.category = "Simple Nucleotide Variation", 
  access = "open", 
  legacy = FALSE, 
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)
GDCdownload(query)
maf <- GDCprepare(query)
maf <- maf %>% read.maf

datatable(getSampleSummary(maf),
          filter = 'top',
          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
          rownames = FALSE)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE)

oncoplot(maf = maf, top = 10, removeNonMutated = TRUE)

# classifies Single Nucleotide Variants into Transitions and Transversions
titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = titv)
