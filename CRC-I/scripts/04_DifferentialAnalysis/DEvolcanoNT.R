# based on Seurat official DE test vignette

library(Seurat)
library(tidyverse)
library(Platypus)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

labelNat <- readRDS("data/nature/CountNat.rds")
labelNat <- ScaleData(labelNat)
labelNat <- SetIdent(labelNat, value = "ITgeno")
# examine its counts and data slots
head(labelNat@assays$RNA@counts)
head(labelNat@assays$RNA@data)
head(labelNat@assays$RNA@scale.data)
# labelNat <- NormalizeData(labelNat)
# IT.markers <- FindMarkers(labelNat, ident.1 = 'IT', fc.name = 'avg_logFC')
# head(IT.markers)
sbclstrNT <- SplitObject(labelNat, split.by = "majorCluster")
# sbclstrNT <- sbclstrNT[c('hB01','hB02','hB03','hB04','hB05','hM02','hM03','hM04','hM12','hM13')]
# 'hB09' cluster has no IT cells
write_rds(sbclstrNT, "data/nature/NatureSubList.rds")

for (i in 1:length(sbclstrNT)) { # 22-23,25-27 have low TT cells
  print(names(sbclstrNT[i]))
  submark <- FindMarkers(sbclstrNT[[i]], ident.1 = "TT", fc.name = "avg_logFC") # use normalized logTPM data
  write.csv(subset(submark, p_val_adj < 0.05), file = paste0("results/DEG/nature/NatureTCount-TTmarkers-", names(sbclstrNT[i]), ".csv"))

  vol <- GEX_volcano2(
    DEGs.input = submark,
    input.type = "findmarkers",
    condition.1 = "TT",
    condition.2 = "II/IT",
    maximum.overlaps = 5
  )

  htmp <- GEX_DEgenes(GEX = sbclstrNT[[i]],
                      min.pct = .25,
                      FindMarkers.out = submark,
                      return.plot = "heatmap",
                      up.genes = 5,
                      down.genes = 5,
                      logFC = F)

  png(filename = paste0("figures/volcano/nature/natureTCount-", names(sbclstrNT[i]), ".png"), width = 600, height = 600)
  plot(vol)
  dev.off()

  png(filename = paste0("figures/heatmap/nature/natureTCount-", names(sbclstrNT[i]), ".png"), width = 600, height = 600)
  plot(htmp[[2]] + theme(text = element_text(size = 16)))
  dev.off()
}

# compare TT and IT only
for (i in 13:length(sbclstrNT)) { # 22,23,25-27 have low TT cells
  print(names(sbclstrNT[i]))
  submark <- FindMarkers(sbclstrNT[[i]],
                         ident.1 = "TT",
                         ident.2 = "IT",
                         fc.name = "avg_logFC") # use normalized logTPM data
  vol <- GEX_volcano2(
    DEGs.input = submark,
    input.type = "findmarkers",
    condition.1 = "TT",
    condition.2 = "IT",
    maximum.overlaps = 5
  )
  png(filename = paste0("figures/volcano/nature/CountIT_TT-", names(sbclstrNT[i]), ".png"), width = 600, height = 600)
  plot(vol)
  dev.off()
  write.csv(subset(submark, p_val_adj < 0.05), file = paste0("results/DEG/nature/CountIT_TTmarkers-", names(sbclstrNT[i]), ".csv"))

  htmp <- GEX_DEgenes(GEX = sbclstrNT[[i]], min.pct = .25, FindMarkers.out = submark, return.plot = "heatmap", group1 = "TT", group2 = "IT", up.genes = 5, down.genes = 5, logFC = F)
  png(filename = paste0("figures/heatmap/nature/natureTT-ITCount-", names(sbclstrNT[i]), ".png"), width = 600, height = 600)
  plot(htmp[[2]] + theme(text = element_text(size = 20)))
  dev.off()
}

GEX_proportions_barplot(GEX = labelNat, stacked.plot = T, source.group = "majorCluster", target.group = "ITgeno")

# auto-annotate cell type by Platypus
labelNat <- GEX_phenotype(labelNat, default = 1)

library(data.table)
meta <- as.data.table(labelNat@meta.data)
meta$UniqueCell_ID
ggplot(meta, aes(x = ITgeno, fill = cell.state)) +
  geom_bar(position = "fill") +
  scale_fill_viridis_d(option = "turbo")

for (i in 14:length(sbclstrNT)) {
  try(submark <- FindMarkers(sbclstrNT[[i]], ident.1 = "TT"))
  genelist <- submark[, 2] # vector of fold change
  names(genelist) <- rownames(submark) # vector of gene id
  genelist <- sort(genelist, decreasing = TRUE)

  gogsea <- gseGO(
    geneList = genelist,
    keyType = "SYMBOL",
    pvalueCutoff = 0.05,
    OrgDb = org.Hs.eg.db
  )
  try(goplot <- dotplot(gogsea))
  png(paste0("figures/GSEAdot/nature/", names(sbclstrNT[i]), ".png"), width = 600, height = 600)
  plot(goplot + theme(text = element_text(size = 16)))
  dev.off()
}

markerset <- read_delim("ref/markerSet.txt")
VlnPlot(labelNat, features = markerset$Exhaustion)
RidgePlot(labelNat, features = markerset$Exhaustion)
RidgePlot(labelNat, features = markerset$Proliferation)
hypoxia <- markerset$Hypoxia
for (i in 1:38) {
  RidgePlot(labelNat, features = hypoxia[i])
  ggsave(paste0("figures/ridgeplot/", hypoxia[i], ".png"))
}
