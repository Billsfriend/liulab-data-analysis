library(tximeta)
library(tidyverse)

coldata <- read_delim("tximeta/coldata.csv", 
                      delim = "\t", escape_double = FALSE, 
                      trim_ws = TRUE)
# SummarizedExperiment object is generated as result
se <- tximeta(coldata)
# transform transcripts-level data to gene level
gse <- summarizeToGene(se)
dim(gse) # counts ~35000 genes here, fewer than 40000+ by hisat2/featureCount workflow
write_rds(gse, 'data/tximetaGSE.rds')

