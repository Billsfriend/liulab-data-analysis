# ref: https://www.jianshu.com/p/adda4536b2cb

library(GSVA)
library(GSVAdata)
library(tidyverse)
library(pheatmap)
library(patchwork)
library(msigdbr)
library(limma) #limma做差异分析

data(c2BroadSets)

m_df = msigdbr(species = "Homo sapiens",
               category = "C2") # C3, C5-GO are worth trying
m_df = subset(m_df, gs_subcat != 'CGP')
# C2 (-CGP) is still too big, just CP:reactome?
msigdbr_list = split(x = m_df$gene_symbol, f = m_df$gs_name)

capCounts <- read_csv('../results/hisat-deg.csv') %>%
  dplyr::select(c('symbol', 'GFP', 'Cap1', 'Cap2'))

expr <- as.matrix(capCounts[2:4])
rownames(expr) <- capCounts$symbol

test <- gsva(expr,
             msigdbr_list,
             kcdf = 'Poisson',
             parallel.sz = 10,
             min.sz = 3)

group <- factor(x = c('GFP', 'Cap1', 'Cap2'),
                levels = c('GFP', 'Cap1', 'Cap2'),
                labels = c('GFP', 'Capsid', 'Capsid'),
                ordered = FALSE)

design <- model.matrix(~ group) # 构建比较矩阵
colnames(design) <- levels(group)
fit <- lmFit(test, design) # fit in linear model
fit1 <- eBayes(fit) # Bayes test
allDiff = topTable(
  fit1,
  adjust = 'fdr',
  coef = 2,
  number = Inf,
  resort.by = 't'
)

allDiff$pathway <- str_to_lower(rownames(allDiff))

write_csv(allDiff, '../results/gsva_results.csv')

# generate heatmap by geom_tile
ggplot(natperc, aes(x=factor,y=subset,fill=mean))+
  geom_tile()+
  geom_raster()+
  scale_fill_viridis_c()+
  labs(fill='Proportion in all T cells')+
  theme(text=element_text(size=16))
