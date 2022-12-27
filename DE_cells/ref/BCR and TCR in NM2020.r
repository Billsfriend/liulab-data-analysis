# edit from Zhang XK's script by Gao J
{
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(cowplot)
  library(pheatmap)
  options(stringsAsFactors = F)
}

# load data
Human_PBMC <- readRDS("covid19/data/blish_covid_NM2020.rds")
DefaultAssay(Human_PBMC) <- "RNA" # set default assay

Idents(Human_PBMC) <- Human_PBMC$Stage
table(Idents(Human_PBMC))

pbmc <-
  subset(x = Human_PBMC, idents = "Ctrl") # only want the control group in this dataset

table(Idents(pbmc)) # make sure only get the control group
head(pbmc, 5)

table(pbmc$cell_type) # make it clear to watch

my_level <- c(
  "B",
  "Class-switched B",
  "IgA PB",
  "IgG PB",
  "CD4 T",
  "CD4m T",
  "CD4n T",
  "CD8eff T",
  "CD8m T",
  "gd T",
  "NK",
  "CD14 Monocyte",
  "CD16 Monocyte",
  "pDC",
  "DC",
  "Neutrophil",
  "SC & Eosinophil",
  "Activated Granulocyte",
  "Platelet",
  "RBC"
)

Idents(pbmc) <- factor(pbmc$cell.type.fine, levels = my_level)
table(Idents(pbmc))

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

subset(pbmc, percent.mt < 10) -> pbmc

pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
DefaultAssay(pbmc) <- "RNA" # set default assay


# BCR --------------
pbmc
{
  features <- read.table("~/ImmR/DE cells/BCR1.txt")
  features <- features$V1
  print(features)
}

features1 <- features[1:10]
features2 <- features[11:16]

features3 <- features[!features %in% c("IGLC6")]

# violin plot
VlnPlot(pbmc, features = features1, ncol = 3) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/NM2020/DE/output/0520/VlnPlot_BCR1.png", width = 15, height = 12)

VlnPlot(pbmc, features = features2, ncol = 3) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/NM2020/DE/output/0520/VlnPlot_BCR2.png", width = 15, height = 12)

VlnPlot(pbmc, features = features3, stack = T, flip = T) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/NM2020/DE/output/0520/VlnPlotstack_BCR1.png", width = 15, height = 12)

# dot plot
DotPlot(pbmc, features = features3) + RotatedAxis() + scale_color_gradientn(colours = rev(brewer.pal(n = 5, name = "RdYlBu")))
ggsave("~/ImmR/NM2020/DE/output/0520/DotPlot_BCR1.png", width = 12, height = 10)

# create subset
{
  select_1 <- subset(x = pbmc, slot = "data", subset = CD79A > 0.0)
  select_2 <- subset(x = pbmc, slot = "data", subset = CD79B > 0.0)
  select_3 <- subset(x = pbmc, slot = "data", subset = IGKC > 0.0)
  select_4 <- subset(x = pbmc, slot = "data", subset = IGLC2 > 0.0)
  select_5 <- subset(x = pbmc, slot = "data", subset = IGLC3 > 0.0)
  select_6 <- subset(x = pbmc, slot = "data", subset = IGLC7 > 0.0)
  select_7 <- subset(x = pbmc, slot = "data", subset = IGHA1 > 0.0)
  select_8 <- subset(x = pbmc, slot = "data", subset = IGHA2 > 0.0)
  select_9 <- subset(x = pbmc, slot = "data", subset = IGHG1 > 0.0)
  select_10 <- subset(x = pbmc, slot = "data", subset = IGHG2 > 0.0)
  select_11 <- subset(x = pbmc, slot = "data", subset = IGHG3 > 0.0)
  select_12 <- subset(x = pbmc, slot = "data", subset = IGHG4 > 0.0)
  select_13 <- subset(x = pbmc, slot = "data", subset = IGHD > 0.0)
  select_14 <- subset(x = pbmc, slot = "data", subset = IGHE > 0.0)
  select_15 <- subset(x = pbmc, slot = "data", subset = IGHM > 0.0)
}

{
  P_1 <- as.data.frame(table(Idents(select_1)))
  P_2 <- as.data.frame(table(Idents(select_2)))
  P_3 <- as.data.frame(table(Idents(select_3)))
  P_4 <- as.data.frame(table(Idents(select_4)))
  P_5 <- as.data.frame(table(Idents(select_5)))
  P_6 <- as.data.frame(table(Idents(select_6)))
  P_7 <- as.data.frame(table(Idents(select_7)))
  P_8 <- as.data.frame(table(Idents(select_8)))
  P_9 <- as.data.frame(table(Idents(select_9)))
  P_10 <- as.data.frame(table(Idents(select_10)))
  P_11 <- as.data.frame(table(Idents(select_11)))
  P_12 <- as.data.frame(table(Idents(select_12)))
  P_13 <- as.data.frame(table(Idents(select_13)))
  P_14 <- as.data.frame(table(Idents(select_14)))
  P_15 <- as.data.frame(table(Idents(select_15)))
}

{
  {
    cell <- as.data.frame(my_level)
    colnames(cell) <- "Var1"

    cell_num <- data.frame(x = table(pbmc@active.ident))
    colnames(cell_num)[1] <- "Var1"
  }

  # total
  test <- join_all(list(
    cell, cell_num, P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9, P_10,
    P_11, P_12, P_13, P_14, P_15
  ), by = "Var1", type = "left")

  test[is.na(test)] <- 0
  colnames(test) <- c("cell", "cell_num", features3)

  test <- test[, 2:length(colnames(test))]

  for (i in 1:nrow(test))
  {
    test[i, ] <- test[i, ] / test[i, 1]
  }

  test <- test[, -1]
  rownames(test) <- my_level

  test <- as.data.frame(t(test))

  # creating a function to compute percentage
  my_percent <- function(num, digits = 2, ...) {
    percentage <- formatC(num * 100, format = "f", digits = digits, ...)

    # appending "%" symbol at the end of
    # calculate percentage value
    paste0(percentage, "%")
  }

  test5 <- as.data.frame(lapply(test, my_percent, 2))
  rownames(test5) <- c(features3)
}

# Generate heat map
hr <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "row", number_color = "black", fontsize = 15,
  labels_row = c(features3), labels_col = my_level
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmaprow_BCR1.png", hr, dpi = 300, width = 15, height = 12)

hc <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "column", number_color = "black", fontsize = 15,
  labels_row = c(features3), labels_col = my_level
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmapcol_BCR1.png", hc, dpi = 300, width = 15, height = 12)

#############################################
# TCR
#############################################
pbmc
{
  features <- read.table("~/ImmR/DE cells/TCR1.txt")
  features <- features$V1
  print(features)
}

features4 <- features[!features %in% c("TRBC1")]

# violin plot
VlnPlot(pbmc, features = features, ncol = 3) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/NM2020/DE/output/0520/VlnPlot_TCR1.png", width = 15, height = 12)

VlnPlot(pbmc, features = features4, stack = T, flip = T) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/NM2020/DE/output/0520/VlnPlotstack_TCR1.png", width = 15, height = 12)

# dot plot
DotPlot(pbmc, features = features4) + RotatedAxis() + scale_color_gradientn(colours = rev(brewer.pal(n = 5, name = "RdYlBu")))
ggsave("~/ImmR/NM2020/DE/output/0520/DotPlot_TCR1.png", width = 12, height = 10)

memory.limit(size = 5000000)

# create subset
{
  select_1 <- subset(x = pbmc, slot = "data", subset = CD3D > 0.0)
  select_2 <- subset(x = pbmc, slot = "data", subset = CD3E > 0.0)
  select_3 <- subset(x = pbmc, slot = "data", subset = CD3G > 0.0)
  select_4 <- subset(x = pbmc, slot = "data", subset = CD247 > 0.0)
  select_5 <- subset(x = pbmc, slot = "data", subset = TRAC > 0.0)
  select_6 <- subset(x = pbmc, slot = "data", subset = TRBC2 > 0.0)
  select_7 <- subset(x = pbmc, slot = "data", subset = TRDC > 0.0)
  select_8 <- subset(x = pbmc, slot = "data", subset = TRGC1 > 0.0)
  select_9 <- subset(x = pbmc, slot = "data", subset = TRGC2 > 0.0)
}

{
  P_1 <- as.data.frame(table(Idents(select_1)))
  P_2 <- as.data.frame(table(Idents(select_2)))
  P_3 <- as.data.frame(table(Idents(select_3)))
  P_4 <- as.data.frame(table(Idents(select_4)))
  P_5 <- as.data.frame(table(Idents(select_5)))
  P_6 <- as.data.frame(table(Idents(select_6)))
  P_7 <- as.data.frame(table(Idents(select_7)))
  P_8 <- as.data.frame(table(Idents(select_8)))
  P_9 <- as.data.frame(table(Idents(select_9)))
}

{
  {
    cell <- as.data.frame(my_level)
    colnames(cell) <- "Var1"

    cell_num <- data.frame(x = table(pbmc@active.ident))
    colnames(cell_num)[1] <- "Var1"
  }

  # total
  test <- join_all(list(cell, cell_num, P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9), by = "Var1", type = "left")

  test[is.na(test)] <- 0
  colnames(test) <- c("cell", "cell_num", features4)

  test <- test[, 2:length(colnames(test))]

  for (i in 1:nrow(test))
  {
    test[i, ] <- test[i, ] / test[i, 1]
  }

  test <- test[, -1]
  rownames(test) <- my_level

  test <- as.data.frame(t(test))

  # creating a function to compute percentage
  my_percent <- function(num, digits = 2, ...) {
    percentage <- formatC(num * 100, format = "f", digits = digits, ...)

    # appending "%" symbol at the end of
    # calculate percentage value
    paste0(percentage, "%")
  }

  test5 <- as.data.frame(lapply(test, my_percent, 2))
  rownames(test5) <- c(features4)
}

# Generate heat map
hr <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "row", number_color = "black", fontsize = 15,
  labels_row = c(features4), labels_col = my_level
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmaprow_TCR1.png", hr, dpi = 300, width = 15, height = 12)

hc <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "column", number_color = "black", fontsize = 15,
  labels_row = c(features4), labels_col = my_level
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmapcol_TCR1.png", hc, dpi = 300, width = 15, height = 12)

#############################################
# DE cells
#############################################

# create subset
pbmc
{
  select_DE_1 <- subset(x = pbmc, slot = "data", subset = CD79A > 0.0 & CD79B > 0.0)
  select_DE_2 <- subset(x = select_DE_1, slot = "data", subset = IGKC > 0.0 | IGLC2 > 0.0 | IGLC3 > 0.0 | IGLC7 > 0.0)
  select_DE_3 <- subset(x = select_DE_2, slot = "data", subset = IGHA1 > 0.0 | IGHA2 > 0.0 | IGHG1 > 0.0 | IGHG2 > 0.0 | IGHG3 > 0.0 | IGHG4 > 0.0 | IGHD > 0.0 | IGHE > 0.0 | IGHM > 0.0)
  select_DE_4 <- subset(x = select_DE_3, slot = "data", subset = CD3D > 0.0 & CD3E > 0.0 & CD3G > 0.0 & CD247 > 0.0)
  select_DE_5 <- subset(x = select_DE_4, slot = "data", subset = (TRAC > 0.0 & (TRBC2 > 0.0)) | (TRDC > 0.0 & (TRGC1 > 0.0 | TRGC2 > 0.0)))
}

{
  P_DE_1 <- as.data.frame(table(Idents(select_DE_1)))
  P_DE_2 <- as.data.frame(table(Idents(select_DE_2)))
  P_DE_3 <- as.data.frame(table(Idents(select_DE_3)))
  P_DE_4 <- as.data.frame(table(Idents(select_DE_4)))
  P_DE_5 <- as.data.frame(table(Idents(select_DE_5)))
}

{
  {
    cell <- as.data.frame(my_level)
    colnames(cell) <- "Var1"

    cell_num <- data.frame(x = table(pbmc@active.ident))
    colnames(cell_num)[1] <- "Var1"
  }

  # total
  test <- join_all(list(cell, cell_num, P_DE_1, P_DE_2, P_DE_3, P_DE_4, P_DE_5), by = "Var1", type = "left")

  test[is.na(test)] <- 0
  colnames(test) <- c("cell", "cell_num", "BCR", "BCR_IGL", "BCR_IGL_IGH", "BCR_IGL_IGH_CD3", "BCR_IGL_IGH_CD3_TCR")

  test <- test[, 2:length(colnames(test))]

  for (i in 1:nrow(test))
  {
    test[i, ] <- test[i, ] / test[i, 1]
  }

  test <- test[, -1]
  rownames(test) <- my_level

  test <- as.data.frame(t(test))

  # creating a function to compute percentage
  my_percent <- function(num, digits = 2, ...) {
    percentage <- formatC(num * 100, format = "f", digits = digits, ...)

    # appending "%" symbol at the end of
    # calculate percentage value
    paste0(percentage, "%")
  }

  test5 <- as.data.frame(lapply(test, my_percent, 2))
  rownames(test5) <- c("BCR", "BCR_IGL", "BCR_IGL_IGH", "BCR_IGL_IGH_CD3", "BCR_IGL_IGH_CD3_TCR")
}

# Generate heat map
hr <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "row", number_color = "black", fontsize = 15,
  labels_row = c("CD79", "CD79_IGL", "CD79_IGL_IGH", "CD79_IGL_IGH_CD3", "CD79_IGL_IGH_CD3_TCR"), labels_col = my_level,
  main = "DE cells (BCR+TCR+) in NM2020 data"
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmaprow_DE_NM2020.pdf", hr, dpi = 300, width = 15, height = 12)


# load data
{
  Human_PBMC <- readRDS("~/blish_covid.seu.rds")
  DefaultAssay(Human_PBMC) <- "RNA" # set default assay

  Idents(Human_PBMC) <- Human_PBMC$Status
  table(Idents(Human_PBMC))

  pbmc <- subset(x = Human_PBMC, idents = "COVID") # only want the control group in this dataset

  table(Idents(pbmc)) # make sure only get the control group
  head(pbmc, 5)

  table(pbmc$cell.type.fine) # make it clear to watch

  my_level <- c(
    "B", "Class-switched B", "IgA PB", "IgG PB", "CD4 T", "CD4m T", "CD4n T", "CD8eff T", "CD8m T", "gd T", "NK", "CD14 Monocyte", "CD16 Monocyte",
    "pDC", "DC", "Neutrophil", "SC & Eosinophil", "Activated Granulocyte", "Platelet", "RBC"
  )

  Idents(pbmc) <- factor(pbmc$cell.type.fine, levels = my_level)
  table(Idents(pbmc))
}

VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
DefaultAssay(pbmc) <- "RNA" # set default assay

#############################################
# DE cells IN COVID
#############################################

# create subset
pbmc
{
  select_DE_1 <- subset(x = pbmc, slot = "data", subset = CD79A > 0.0 & CD79B > 0.0)
  select_DE_2 <- subset(x = select_DE_1, slot = "data", subset = IGKC > 0.0 | IGLC2 > 0.0 | IGLC3 > 0.0 | IGLC7 > 0.0)
  select_DE_3 <- subset(x = select_DE_2, slot = "data", subset = IGHA1 > 0.0 | IGHA2 > 0.0 | IGHG1 > 0.0 | IGHG2 > 0.0 | IGHG3 > 0.0 | IGHG4 > 0.0 | IGHD > 0.0 | IGHE > 0.0 | IGHM > 0.0)
  select_DE_4 <- subset(x = select_DE_3, slot = "data", subset = CD3D > 0.0 & CD3E > 0.0 & CD3G > 0.0 & CD247 > 0.0)
  select_DE_5 <- subset(x = select_DE_4, slot = "data", subset = (TRAC > 0.0 & (TRBC2 > 0.0)) | (TRDC > 0.0 & (TRGC1 > 0.0 | TRGC2 > 0.0)))
}

{
  P_DE_1 <- as.data.frame(table(Idents(select_DE_1)))
  P_DE_2 <- as.data.frame(table(Idents(select_DE_2)))
  P_DE_3 <- as.data.frame(table(Idents(select_DE_3)))
  P_DE_4 <- as.data.frame(table(Idents(select_DE_4)))
  P_DE_5 <- as.data.frame(table(Idents(select_DE_5)))
}

{
  {
    cell <- as.data.frame(my_level)
    colnames(cell) <- "Var1"

    cell_num <- data.frame(x = table(pbmc@active.ident))
    colnames(cell_num)[1] <- "Var1"
  }

  # total
  test <- join_all(list(cell, cell_num, P_DE_1, P_DE_2, P_DE_3, P_DE_4, P_DE_5), by = "Var1", type = "left")

  test[is.na(test)] <- 0
  colnames(test) <- c("cell", "cell_num", "BCR", "BCR_IGL", "BCR_IGL_IGH", "BCR_IGL_IGH_CD3", "BCR_IGL_IGH_CD3_TCR")

  test <- test[, 2:length(colnames(test))]

  for (i in 1:nrow(test))
  {
    test[i, ] <- test[i, ] / test[i, 1]
  }

  test <- test[, -1]
  rownames(test) <- my_level

  test <- as.data.frame(t(test))

  # creating a function to compute percentage
  my_percent <- function(num, digits = 2, ...) {
    percentage <- formatC(num * 100, format = "f", digits = digits, ...)

    # appending "%" symbol at the end of
    # calculate percentage value
    paste0(percentage, "%")
  }

  test5 <- as.data.frame(lapply(test, my_percent, 2))
  rownames(test5) <- c("BCR", "BCR_IGL", "BCR_IGL_IGH", "BCR_IGL_IGH_CD3", "BCR_IGL_IGH_CD3_TCR")
}

# Generate heat map
hr <- pheatmap(as.matrix(test),
  cluster_row = FALSE, cluster_cols = FALSE,
  display_numbers = test5,
  scale = "row", number_color = "black", fontsize = 15,
  labels_row = c("CD79", "CD79_IGL", "CD79_IGL_IGH", "CD79_IGL_IGH_CD3", "CD79_IGL_IGH_CD3_TCR"), labels_col = my_level,
  main = "DE cells (BCR+TCR+) in COVID from NM2020 data"
)
ggsave("~/ImmR/NM2020/DE/output/0520/heatmaprow_DE_NM2020_COVID.pdf", hr, dpi = 300, width = 15, height = 12)
