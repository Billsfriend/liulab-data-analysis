library(edgeR)
library(tidyverse)
library(data.table)
library(clusterProfiler)
# edgeRUsersGuide()

hisatIgM <- read_delim('data/Galaxy11-[featureCounts_IgM__Counts].tabular')
hisatIgG1 <- read_delim('data/Galaxy19-[featureCounts_IgG__Counts].tabular')

hisatIgM$IgG1 <- hisatIgG1$counts

counts <- data.frame(hisatIgM[,-1])

rownames(counts) <- hisatIgM$Geneid

edglist <- DGEList(counts = counts, group = factor(c('IgM','IgG1')))

# filter out low expression genes
keep <- filterByExpr(edglist)
edglist <- edglist[keep, , keep.lib.sizes=FALSE]

# normalize for effective library size
edglist <- calcNormFactors(edglist)

# estimate dispersion
# edglist <- estimateDisp(edglist)

# for samples without replicates, we may choose arbitrary BCV (square-root-dispersion) value
# typically, in well-controlled experiment, the value is 0.4 for human specimen, 0.1 for genetically identical model organisms, and 0.01 for technical replicates
bcv <- 0.15

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

symbols<- bitr(rownames(allDeg), fromType = 'ENTREZID', toType = 'SYMBOL', OrgDb = 'org.Mm.eg.db', drop = F)

allDeg$symbol <- symbols$SYMBOL
allDeg$entrez <- symbols$ENTREZID

# manually calculate adjusted P value
allDeg$padj <- p.adjust(allDeg$PValue, method = 'BH')

# calculate count-per-million for both group
cpm <- as_tibble(cpmByGroup(edglist))

allDeg$IgM <- cpm$IgM
allDeg$IgG1 <- cpm$IgG1

fwrite(allDeg, 'results/hisat-deg-entrez.csv')

# plotMD logFC-AveLogCPM
plotMD(LRTreat, main = 'IgG-B cell differential expressed genes compared with IgM-B cell')

# gene set testing
# fry(edglist, dispersion = bcv^2)
# not allowed in no-replicates...

# even alternative splicing analysis is possible
# diffSpliceDGE(fit)
