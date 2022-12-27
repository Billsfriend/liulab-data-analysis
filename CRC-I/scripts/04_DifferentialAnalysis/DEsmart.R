# based on Seurat official DE test vignette

library(Seurat)
library(tidyverse)
library(Platypus)

labelSmart <- readRDS("data/cell/labelSmart.rds")

labelSmart <- ScaleData(labelSmart)
labelSmart <- SetIdent(labelSmart, value = 'ITgeno')
#IT.markers <- FindMarkers(labelSmart, ident.1 = 'IT', fc.name = 'avg_logFC')
#head(IT.markers)
subclstrSmart <- SplitObject(labelSmart, split.by = 'Sub_Cluster')
# subclstrSmart <- subclstrSmart[c('hB01','hB02','hB03','hB04','hB05','hM02','hM03','hM04','hM12','hM13')]

write_rds(subclstrSmart, 'data/cell/subclstrSmart.rds')
#subclstrSmart <- readRDS('data/cell/subclstrSmart.rds')

for (i in 29:length(subclstrSmart)) { # cluster 28,44-48 has low IT cells
  cname = names(subclstrSmart[i])
  print(cname)
  submark <- FindMarkers(subclstrSmart[[i]], ident.1 = 'IT', fc.name = 'avg_logFC')
  vol <- GEX_volcano(DEGs.input = submark,
              input.type = 'findmarkers',
              condition.1 = 'IT',
              condition.2 = 'II', maximum.overlaps = 15)
  png(filename = paste0('figures/volcano/smart/Smart-', cname, '.png'), width = 600, height = 600)
  plot(vol)
  dev.off()
  write.csv(subset(submark, p_val_adj < 0.05), file = paste0('results/DEG/smart/cell-Smart-ITmarkers-', cname, '.csv'))
}


