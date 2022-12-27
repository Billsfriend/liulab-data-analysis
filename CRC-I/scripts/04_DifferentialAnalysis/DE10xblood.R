# based on Seurat official DE test vignette

library(Seurat)
library(tidyverse)
library(Platypus)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
source('Scripts/functions/GEX_volcano.R')
source('Scripts/functions/DEG.R')

# label10x <- readRDS("data/cell/label10x.rds")
# subclstr10x <- SplitObject(label10x, split.by = 'Tissue')
# write_rds(subclstr10x[['P']], 'data/cell/blood10x.rds')
# write_rds(subclstr10x[['N']], 'data/cell/adjcent10x.rds')

label10x <- read_rds('data/cell/blood10x.rds')
subclstr10x <- SplitObject(label10x, split.by = 'Sub_Cluster')

# prepare to get DEGs numbers for each subcluster
#degnum <- tibble(id=1:length(subclstr10x), cellset = NA, degs = NA)

for (i in 1:length(subclstr10x)) {
  
  res <- tryCatch(deg(subclstr10x[i], id1='IT',id2='II'),
           error = function(e)
             return(NULL))
  if (is.null(res)){
    next
  }
  print(res[[1]])
  print(paste0(i, '/', length(subclstr10x)))
 
   write.csv(subset(res[[2]], p_val_adj < 0.05), file = paste0('results/DEG/blood10x/10xTumor-ITmarkers-',res[[1]], '.csv'))
  write.csv(subset(res[[2]]), file = paste0('results/DEG/blood10x/ALL10xTumor-ITmarkers-',res[[1]], '.csv'))
  png(filename = paste0('figures/volcano/blood10x/10xTumor-', res[[1]], '.png'), width = 600, height = 600)
  plot(res[[3]] +
         theme(text = element_text(size = 14))
  )
  dev.off()
  
  png(filename = paste0('figures/heatmap/10x/10xBlood-IT-', res[[1]], '.png'), width = 600, height = 600)
  plot(res[[4]]+
         scale_fill_viridis_c()+
         theme(text = element_text(size = 20)))
  dev.off()
  
  data <- res[[2]]
  genelist <- data[,2] # vector of fold change
  names(genelist) = rownames(data)# vector of gene id
  genelist = sort(genelist, decreasing = TRUE)
  
  gogsea <- gseGO(
    geneList = genelist,
    keyType = 'SYMBOL',
    eps=0,
    OrgDb = org.Hs.eg.db
  )
  rs <- tryCatch(dotplot(gogsea, showCategory=10),
           error=function(e)
             return(NULL))
  if(is.null(rs)){
    next
  }
  png(filename = paste0('figures/GSEAdot/10x/10xBlood-IT-', res[[1]], '.png'), width = 600, height = 600)
  plot(dotplot(gogsea, showCategory=10) + ggtitle(paste0(res[[1]], ' IT Gene Ontology')))
  dev.off()
  
}

# bar plot: number of DEGs from each subcluster pairs
degnum <- subset(degnum, !is.na(degs) & degs != 0) # get rid of pairs with no significant DEGs

submark <- FindMarkers(nk16, ident.1 = 'IT', fc.name = 'avg_logFC', logfc.threshold = 0.1)

vol <- GEX_volcano(DEGs.input = submark,
                    input.type = 'findmarkers',
                    condition.1 = 'IT',
                    n.label.down = 6,
                    n.label.up = 6,
                    maximum.overlaps = 6)
plot(vol+
       theme(text = element_text(size = 14)))

unilist <- rownames(nk16@assays$RNA@counts)


nkup <- subset(submark, p_val_adj < 0.05)
genelist <- nkup[,2]
names(genelist) <- rownames(nkup)

oraGo <- enrichGO(gene = names(genelist),
                  keyType = 'SYMBOL',
                  OrgDb = org.Hs.eg.db,
                  universe = unilist,
                  ont = 'ALL',
                  minGSSize = 10,
                  readable = TRUE)

dotplot(oraGo)

heatplot(oraGo, showCategory = 10, foldChange = genelist)
# why no foldchange color?

oraGot <- pairwise_termsim(oraGo) # get similarity matrix
treeplot(oraGot)

VlnPlot(nk16, features = c('PRF1','LCK','NCR3','FCGR3A'))
RidgePlot(nk16, features = c('PRF1','LCK','NCR3','FCGR3A'))
FeaturePlot(nk16, features = c('PRF1','FCGR3A'), split.by = 'ITgeno')
# define major clusters
degnum$bigset <- NA
degnum$bigset <- str_trunc(degnum$cellset, 2, 'right', ellipsis = '')

setDT(degnum)

ggplot(degnum)+
  geom_col(aes(x=cellset, y=degs, fill=bigset))+
  coord_flip()+
  ylab('Number of significant DEGs between II and IT')+
  scale_fill_discrete(
    name = 'Major clusters',
    limits = c("hT", "hM", "hI", 'hB'),
    labels = c("T cell", "Myeloid cell", "ILC", 'B cell')
  )

ggplot(degnum[bigset=='hT'])+
  geom_col(aes(x=cellset, y=degs, fill=bigset))+
  coord_flip()+
  ylab('Number of significant DEGs between II and IT')+
  scale_fill_discrete(
    name = 'Major clusters',
    labels = c('T cell')
  )+
  theme(text = element_text(size=16))
