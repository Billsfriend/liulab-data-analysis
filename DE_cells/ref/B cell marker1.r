{
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(magrittr)
  library(ggplot2)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(RColorBrewer)
  library(cowplot)
  library(plyr)
  library(pheatmap)
  options(stringsAsFactors = F)
}

#load data
{
  Human_PBMC <- readRDS("~/Final_nCoV_0716_upload.RDS")
  DefaultAssay(Human_PBMC) <- "RNA"               #set default assay
  pbmc <- subset(x = Human_PBMC, idents = "Ctrl") #only want the control group in this dataset
  
  my_level <- c("Naive B cells","Memory B cells","Plasma","Cycling Plasma","Naive T cells","Activated CD4 T cells","Cytotoxic CD8 T cells","Cycling T cells","MAIT","NKs","XCL+ NKs","Monocytes","Megakaryocytes","DCs","Stem cells")
  
  table(pbmc$cell_type)           #make it clear to watch
  
  Idents(pbmc) <- factor(pbmc$cell_type, levels= my_level)
  
  table(Idents(pbmc))             #make sure only get the control group
}

#############################################
#B cell marker
#############################################
{
  B1 <- read.table("~/ImmR/DE cells/B cell marker1.txt")
  features <- B1$V1
  print(features)
}

#violin plot
VlnPlot(pbmc, features = features) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/DE cells/output/0225_B_cell_marker/VlnPlot1.png",width = 15,height = 12)

VlnPlot(pbmc, features = features, stack = T, flip = T) & theme(axis.title.x = element_blank())
ggsave("~/ImmR/DE cells/output/0225_B_cell_marker/VlnPlotstack1.png",width = 15,height = 12)


#dot plot
DotPlot(pbmc, features = features) + RotatedAxis() + scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
ggsave("~/ImmR/DE cells/output/0225_B_cell_marker/DotPlot.png",width = 12,height = 10)

#UMAP plot
i=1
while (i <= length(features))
{
  a = i
  b = i + 1
  features1 <- features[a:b]
  print(features1)
  FeaturePlot(pbmc, features = features1, 
              label=pbmc$cell_type, 
              cols=c("lightgrey","#DE1F1F"), 
              label.size = 5, 
              pt.size = 0.5)
  ggsave(paste0("~/ImmR/DE cells/output/0225_B_cell_marker/","FeaturePlot",{i},".png"),width = 20,height = 10)
  i = i + 2
  if(a == length(features))
  {
    FeaturePlot(pbmc, features = features[[a]], 
                label=pbmc$cell_type, 
                cols=c("lightgrey","#DE1F1F"), 
                label.size = 5, 
                pt.size = 0.5)
    ggsave(paste0("~/ImmR/DE cells/output/0225_B_cell_marker/","FeaturePlot",{a},".png"),width = 10,height = 10)
  }
}

#create subset
{
  select_1 <- subset(x = pbmc, slot="data",subset = BCL2 > 0.0)
  select_2 <- subset(x = pbmc, slot="data",subset = BCL6 > 0.0)
  select_3 <- subset(x = pbmc, slot="data",subset = CD1D > 0.0)
  select_4 <- subset(x = pbmc, slot="data",subset = CD5 > 0.0)
  select_5 <- subset(x = pbmc, slot="data",subset = CD19 > 0.0)
  select_6 <- subset(x = pbmc, slot="data",subset = CD22 > 0.0)
  select_7 <- subset(x = pbmc, slot="data",subset = CD24 > 0.0)
  select_8 <- subset(x = pbmc, slot="data",subset = CD27 > 0.0)
  select_9 <- subset(x = pbmc, slot="data",subset = CD34 > 0.0)
  select_10 <- subset(x = pbmc, slot="data",subset = CD38 > 0.0)
  select_11 <- subset(x = pbmc, slot="data",subset = CD40 > 0.0)
  select_12 <- subset(x = pbmc, slot="data",subset = CD40LG > 0.0)
  select_13 <- subset(x = pbmc, slot="data",subset = CD53 > 0.0)
  select_14 <- subset(x = pbmc, slot="data",subset = CD69 > 0.0)
  select_15 <- subset(x = pbmc, slot="data",subset = CD72 > 0.0)
  select_16 <- subset(x = pbmc, slot="data",subset = CD80 > 0.0)
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
  P_16 <- as.data.frame(table(Idents(select_16)))
}

{
  {
    cell <- as.data.frame(my_level)
    colnames(cell) <- "Var1"
    
    cell_num <- data.frame(x=table(pbmc$cell_type))
    colnames(cell_num)[1] <- "Var1"
  }
  
  #total
  test <- join_all(list(cell, cell_num, P_1, P_2, P_3, P_4, P_5, P_6, P_7, P_8, P_9, P_10,
                        P_11, P_12, P_13, P_14, P_15, P_16), by='Var1', type='left')
  
  test[is.na(test)] <- 0
  colnames(test) <- c("cell", "cell_num", features)
  
  test <- test[,2:length(colnames(test))]
  
  for (i in 1:nrow(test))
  {
    test[i,] <- test[i,]/test[i,1]
  }
  
  test <- test[,-1]
  rownames(test) <- my_level
  
  test <- as.data.frame(t(test))
  
  # creating a function to compute percentage
  my_percent <- function(num, digits = 2, ...) 
  {      
    percentage <-formatC(num * 100, format = "f", digits = digits, ...)
    
    # appending "%" symbol at the end of
    # calculate percentage value
    paste0(percentage, "%")
  }
  
  test5 <- as.data.frame(lapply(test,my_percent,2))
  rownames(test5) <- c(features)
  
}

#Generate heat map 
hr = pheatmap(as.matrix(test),
              cluster_row = FALSE,cluster_cols = FALSE,
              display_numbers = test5,
              scale="row",number_color = "black",fontsize = 18, 
              labels_row =c(features), labels_col = my_level)
ggsave("~/ImmR/DE cells/output/0225_B_cell_marker/heatmaprow1.png", hr, dpi = 300, width = 15, height = 12)

hc = pheatmap(as.matrix(test),
              cluster_row = FALSE,cluster_cols = FALSE,
              display_numbers = test5,
              scale="column",number_color = "black",fontsize = 18, 
              labels_row =c(features), labels_col = my_level)
ggsave("~/ImmR/DE cells/output/0225_B_cell_marker/heatmapcol1.png", hc, dpi = 300, width = 15, height = 12)

