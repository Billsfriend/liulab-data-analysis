library(edgeR)
library(tidyverse)
library(data.table)

edglist <- read_rds('data/edglist.rds')

# filter out low expression genes
keep <- filterByExpr(edglist, group = 1:2)
edglist <- edglist[keep, , keep.lib.sizes=FALSE]

# normalize for effective library size
# edglist <- calcNormFactors(edglist) # offsets pre-calculated by tximeta, need not to re-compute

# estimate dispersion
# edglist <- estimateDisp(edglist)

# for samples without replicates, we may choose arbitrary BCV (square-root-dispersion) value
# typically, in well-controlled experiment, the value is 0.4 for human specimen, 0.1 for genetically identical model organisms, and 0.01 for technical replicates
bcv <- 0.15

edglist$samples$group <- factor(c('IgM','IgG'))

# perform exact test
et15 <- exactTest(edglist, dispersion = bcv^2)

SDEG <- topTags(et15)

# or, fit in GLM (generalized linear model)
fit <- glmFit(edglist, dispersion = bcv^2)
# and perform likelihood-ratio test
LRT <- glmLRT(fit)

SDEG.lrt <- topTags(LRT)
# for bulk-seq with more replicates, quasi-likelihood F-test is more recommended: glmQLFit(), glmQLFTest()

# edgeR provide an improved TREAT test for a better threshold to detect SDEGs
LRTreat <- glmTreat(fit)

# summarized up/down-regulated SDEG numbers
testRes <- decideTests(LRTreat)
summary(testRes)

allDeg <- LRTreat$table
allDeg$symbol <- LRTreat$genes$gene_name
allDeg$ensembl <- rownames(allDeg) # easier to read
allDeg$entrez <- LRTreat$genes$entrezid # use in clusterProfiler

# manually calculate adjusted P value
allDeg$padj <- p.adjust(allDeg$PValue, method = 'BH')

# calculate count-per-million for both group
cpm <- as_tibble(cpmByGroup(edglist))

allDeg$IgM <- cpm$IgM
allDeg$IgG <- cpm$IgG

# extract 700 most significant degs
SDEG.lrTreat <- topTags(LRTreat, n=700)

sigDeg <- SDEG.lrTreat$table

write_csv(sigDeg, 'results/sig-deg.csv')

fwrite(allDeg, 'results/all-deg-entrez.csv')

# plotMD logFC-AveLogCPM
plotMD(LRTreat, main = 'IgG-B cell differential expressed genes compared with IgM-B cell')

# gene set testing
# fry(edglist, dispersion = bcv^2)
# not allowed in no-replicates...

# even alternative splicing analysis is possible
# diffSpliceDGE(fit)
