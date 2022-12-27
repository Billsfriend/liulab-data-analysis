library(GSVA)
library(GSVAdata)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(msigdbr)
library(limma) #limma做差异分析

data(c2BroadSets)

# read DESeq normalized results
data <- read.csv("DiffAnalysis/DEG_SalTx_results.all.csv", row.names = 1)
# only keep normalized expression data and make it matrix
normalExpr <- as.matrix(subset(data, select = c(A1,A2,A3,B1,B2,B3)))
logNorExpr <- log10(normalExpr+1)
# also read gene counts data for Poisson model calc
counts <- as.matrix(read.csv('data/SalTxiCounts.csv', row.names = 1))

# level design
meta <- c('con', 'con', 'con', 'cko', 'cko', 'cko')
group <-
  factor(meta, levels = c("con", "cko"), ordered = F)## put control group first
design <- model.matrix( ~ group) # construct a matrix for comparison
colnames(design) <- levels(group)

# select gene set list of interest
m_df = msigdbr(species = "Mus musculus", category = "C3")
m_df = subset(m_df, gs_subcat != 'CGP')
# H Hallmark
# C2 Curated genes
# C3 Regulatory targets
# C5 Ontology genes
msigdbr_list = split(x = m_df$ensembl_gene, f = m_df$gs_name) # expr use ensembl name as row names

# perform gene set variation analysis
CanonP <- gsva(logNorExpr, msigdbr_list, kcdf='Gaussian', parallel.sz = 10, min.sz = 3)

fit <- lmFit(CanonP, design)#线性模型拟合
fit1 <- eBayes(fit)#贝叶斯检验
allDiff = topTable(
  fit1,
  adjust = 'fdr',
  coef = 2,
  number = Inf,
  resort.by = 't',
  p.value = 0.05
)
write.csv(allDiff, "results/GSVA-ReguTaget-ST.csv")
