library(edgeR)
library(clusterProfiler)
library(tidyverse)
library(pheatmap)
library(org.Hs.eg.db)
library(PoiClaClu)
# edgeRUsersGuide()

# import data tidy
countGFP <- read_delim('../data/RNAseq/GFP_Counts.tabular')
countCap1 <- read_delim('../data/RNAseq/Cap1_Counts.tabular')
countCap2 <- read_delim('../data/RNAseq/Cap2_Counts.tabular')

countGFP %>%
  left_join(countCap1) %>%
  left_join(countCap2) ->
  countsData

colnames(countsData) <- c('entrez', 'GFP', 'Cap1', 'Cap2') 

symbols <- bitr(countsData$entrez,
               fromType = 'ENTREZID',
               toType = 'SYMBOL',
               OrgDb = 'org.Hs.eg.db',
               drop = F)

countsData$symbol <- symbols$SYMBOL

write_csv(countsData, '../data/rna_seq_counts.csv')

# create DGEList object --------------
counts <- data.frame(countsData[c(2,4)])

rownames(counts) <- countsData$entrez

edglist <- DGEList(counts = counts,
                   group = c('GFP', 'Capsid'))

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

# for [1 vs 2] group, it is ok to evaluate dispersion
# fit in GLM (generalized linear model)
edglist %>%
  glmFit(dispersion = bcv^2) %>%
  glmLRT() %>% # and perform likelihood-ratio test
  topTags() ->
  glm_lrt_tags

# for bulk-seq with more replicates, quasi-likelihood F-test is more recommended: glmQLFit(), glmQLFTest()

# edgeR provide an improved TREAT test for a better threshold to detect SDEGs
edglist %>%
  glmFit(dispersion = bcv^2) %>%
  glmTreat() ->
  glm_treat_capsid

# summarized up/down-regulated SDEG numbers
testRes <- decideTests(glm_treat_capsid)
summary(testRes)

allDeg <- as_tibble(glm_treat_capsid$table,
                    rownames = 'entrez')

symbols <- bitr(allDeg$entrez,
                fromType = 'ENTREZID',
                toType = 'SYMBOL',
                OrgDb = 'org.Hs.eg.db',
                drop = FALSE)

allDeg$symbol <- symbols$SYMBOL

# manually calculate adjusted P value
allDeg$padj <- p.adjust(allDeg$PValue,
                        method = 'BH')

# calculate count-per-million for both group
edglist %>%
  cpmByGroup(dispersion = bcv^2) %>%
  as_tibble() ->
  cpm

allDeg$GFP_cpm <- cpm$GFP
allDeg$Capsid_cpm <- cpm$Capsid

allDeg$logFC <- -allDeg$logFC

write_csv(allDeg, '../results/cap2-deg.csv')

# visualization -----------
# plotMD logFC-AveLogCPM
plotMD(glm_treat_capsid,
       main = 'Differential expressed genes in 293T cells transfected by GFP or Capsid')

# highlight significant DEGs
allDeg %>%
  filter(padj < 0.01 & abs(logFC) > 1) %>%
  mutate(type = logFC > 0) %>%
  right_join(allDeg) ->
  allDeg

# volcano plot
allDeg %>%
  ggplot(aes(x = logFC, y = -log10(PValue))) +
  geom_point(aes(color = type)) +
  geom_hline(aes(yintercept = 3.9), linetype = "dashed") +
  geom_vline(aes(xintercept = -1), linetype = "dashed") +
  geom_vline(aes(xintercept = 1), linetype = "dashed") +
  scale_color_discrete(labels = c('Enrich in GFP', 'Enrich in Capsid', 'NS'),
                       type = c('blue', 'red', 'grey'))

# heatmap
allDeg %>%
  filter(!is.na(type)) ->
  sigDeg

prettyMat <- as.matrix(sigDeg[11:13])
rownames(prettyMat) <- sigDeg$symbol

pheatmap(prettyMat,
         scale = 'row',
         show_rownames = FALSE,
         angle_col = 0)

# heatmap of IFNG pathway
IFNcounts <- countsData %>%
  filter(symbol %in% ifngPathway) %>%
  filter(GFP + Cap1 + Cap2 > 0) %>%
  filter(GFP * 2 > (Cap1 + Cap2)*2)

ifnMat <- as.matrix(IFNcounts[2:4])
rownames(ifnMat) <- IFNcounts$symbol

pheatmap(ifnMat,
         scale = 'row',
         #show_rownames = FALSE,
         angle_col = 0)

# poisson distance matrix of samples
t <- t(countsData[2:4])
PoissonDis <- PoissonDistance(t)
PoissonDis$dd
