# this script is meant to work with snakemake
# BiocManager::install('tximeta')
library(tximeta)
library(org.Mm.eg.db)

# provide paths of quant.sf and corresponding names
coldata <- read.csv('quants/coldata.csv')
coldata$condition <- as.factor(coldata$condition)

# SummarizedExperiment object is generated as result
se <- tximeta(coldata)

# transform transcripts-level data to gene level
gse <- summarizeToGene(se)

dim(gse) # counts ~35000 genes here, fewer than 40000+ by hisat2/featureCount workflow

# directly make DGEList object for downstream edgeR analysis
edglist <- makeDGEList(gse)

saveRDS(gse, snakemake@output[[1]])
saveRDS(edglist, 'data/edglist.rds')
