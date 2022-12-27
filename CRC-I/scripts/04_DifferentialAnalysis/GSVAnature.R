# ref: https://www.jianshu.com/p/adda4536b2cb

library(GSVA)
library(GSVAdata)
library(Seurat)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(msigdbr)
library(limma) #limma做差异分析

data(c2BroadSets)

sbclstrNT <- readRDS("data/nature/NatureSubList.rds")

m_df = msigdbr(species = "Homo sapiens", category = "C2") # C3, C5-GO are worth trying
m_df = subset(m_df, gs_subcat != 'CGP')
# m_df2 = msigdbr(category = "C5")
# m_df2 = subset(m_df2, gs_subcat != 'HPO')
# m_df = merge(m_df, m_df2, all = 1)
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

for (i in 24:length(sbclstrNT)){ # 23,25,26 have low TT. 27 has low IT
  scl <- sbclstrNT[[i]]
  cname <- names(sbclstrNT[i])
  print(cname)
  # compare 2 groups
  DD1 <-  subset(scl, idents = c("IT", "TT")) # control group first
  expr <- as.data.frame(DD1@assays$RNA@counts)
  meta <- DD1@meta.data[, c("ITgeno")]
  expr = as.matrix(expr)
  kegg <- gsva(expr, msigdbr_list,
               parallel.sz = 10, min.sz = 10, max.sz = 100)
  group <-
    factor(meta, levels = c("IT", "TT"), ordered = F)## put control group first
  design <- model.matrix( ~ group)# 构建比较矩阵
  colnames(design) <- levels(group)
  fit <- lmFit(kegg, design)#线性模型拟合
  fit1 <- eBayes(fit)#贝叶斯检验
  allDiff = topTable(
    fit1,
    adjust = 'fdr',
    coef = 2,
    number = Inf,
    resort.by = 't',
    p.value = 0.05
  )
  write.csv(allDiff, paste0("results/GSVA/Nature-", cname, '.csv'))
  # scripts below are meant to visualize the results
  # up <- head(rownames(allDiff), n=3)
  # down <- tail(rownames(allDiff), n=3)
  # TEST <- c(up, down)
  # p <- allDiff
  # p$ID <- rownames(p)
  # q <- p[TEST, ]
  # group1 <- c(rep("TT", 3), rep("IT", 3))
  # df <- data.frame(ID = q$ID,
  #                  score = q$t,
  #                  group = group1)
  # # 按照score排序
  # sortdf <- df[order(df$score), ]
  # sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列
  # p <- ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') +
  #   coord_flip() +
  #   theme_bw() + #去除背景色
  #   theme(panel.grid = element_blank()) +
  #   theme(panel.border = element_rect(size = 0.6))
  # png(paste0('figures/GSVA/nature-', cname, '.png'),
  #     height = 480,
  #     width = 480)
  # plot(p)
  # dev.off()
}

# scripts below only generate strict GSVA csv
for (i in 24:length(sbclstrNT)){ #23,25,26 have low TT. 27 has low IT
  scl <- sbclstrNT[[i]]
  cname <- names(sbclstrNT[i])
  print(cname)
  # compare 2 groups
  DD1 <-  subset(scl, idents = c("IT", "TT")) # control group first
  expr <- as.data.frame(DD1@assays$RNA@data) # use logTPM data
  meta <- DD1@meta.data[, c("ITgeno")]
  expr = as.matrix(expr)
  kegg <- gsva(expr, msigdbr_list,
               parallel.sz = 10, min.sz = 10)
  group <-
    factor(meta, levels = c("IT", "TT"), ordered = F)## put control group first
  design <- model.matrix( ~ group)# 构建比较矩阵
  colnames(design) <- levels(group)
  fit <- lmFit(kegg, design)#线性模型拟合
  fit1 <- eBayes(fit)#贝叶斯检验
  allDiff = topTable(
    fit1,
    adjust = 'fdr',
    coef = 2,
    number = Inf,
    resort.by = 't',
    p.value = 0.05
  )
  write.csv(allDiff, paste0("results/GSVA/CP/StctNature-", cname, '.csv'))
}
