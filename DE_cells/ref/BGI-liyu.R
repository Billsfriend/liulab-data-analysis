getwd()
{
  library(Seurat)
  library(dplyr)
  library(plyr)
  library(Matrix)
  library(magrittr)
  library(ggplot2)
  library(ggpubr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(RColorBrewer)
  library(cowplot)
  library(pheatmap)
  options(stringsAsFactors = F)
}


#20210919
Human_PBMC <- readRDS("./data/Final_nCoV_0716_upload.RDS")
DefaultAssay(Human_PBMC) <- "RNA"
pbmc<-subset(x = Human_PBMC, idents = "Ctrl")

table(pbmc$cell_type)
table(Idents(pbmc))

all.genes <- rownames(pbmc@assays$RNA@data)
my_level <- c("Naive B cells",
              "Memory B cells",
              "Plasma",
              "Cycling Plasma",
              "Naive T cells",
              "Activated CD4 T cells",
              "Cytotoxic CD8 T cells",
              "Cycling T cells",
              "MAIT",
              "NKs",
              "XCL+ NKs",
              "Monocytes",
              "Megakaryocytes",
              "DCs",
              "Stem cells")
Idents(pbmc) <- factor(pbmc$cell_type, levels = my_level)


#Plot 1
features <- c("CD19","CR2","CD81")
#features <- c("CD4","CD8A")
FeaturePlot(pbmc,
            features = features,
            label=pbmc$cell.type,
            cols=c("lightgrey","#DE1F1F"),
            label.size = 4,
            pt.size = 0.2)
ggsave("./output/1011/FeaturePlot.png",
       width = 10,
       height = 10)
#DotPlot(pbmc, features = features) + RotatedAxis()
DotPlot(pbmc, features = features) +
  RotatedAxis()+
  scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
ggsave("./output/1011/DotPlot.png",
       width = 6,
       height = 8)
DotPlot(pbmc, features = features) +
  RotatedAxis()
ggsave("./output/1011/DotPlot2.png",
       width = 6,
       height = 8)
VlnPlot(object = pbmc, features = features)
ggsave("./output/1011/VlnPlot.png",
       width = 12,
       height = 10)


# CD19 CR2 CD81 table
select_CD19 <- subset(x = pbmc, slot="data",subset = CD19 > 0.0)
select_CR2 <- subset(x = pbmc, slot="data",subset = CR2 > 0.0)
select_CD81 <- subset(x = pbmc, slot="data",subset = CD81 > 0.0)
select_all <- subset(x = pbmc, slot="data",subset = CD19 > 0.0 & CR2 > 0.0 & CD81 > 0.0)

select_PTPRC <- subset(x = pbmc, slot="data",subset = PTPRC > 0.0)
PTPRC_num <- data.frame(table(Idents(select_PTPRC)))
cell_num <- data.frame(x=table(pbmc$cell_type))
colnames(cell_num)[1] <- "Var1"

cell_count <- data.frame()
P_19 <- as.data.frame(table(Idents(select_CD19)))
P_21 <- as.data.frame(table(Idents(select_CR2)))
P_81 <- as.data.frame(table(Idents(select_CD81)))
P_all <- as.data.frame(table(Idents(select_all)))

cell <- as.data.frame(my_level)
colnames(cell) <- "Var1"
test <- left_join(cell,cell_num,by="Var1")
test <- left_join(test,P_19,by="Var1")
test <- left_join(test,P_21,by="Var1")
test <- left_join(test,P_81,by="Var1")
test <- left_join(test,P_all,by="Var1")
test <- left_join(test,PTPRC_num,by="Var1")

colnames(test) <- c("cell","cell_num","CD19","CD21","CD81","CD19+CD21+CD81","CD45")
write.csv(x = test,file = paste0("./output/1011//CD19_CD21_CD81_CD45count.csv"))
test[is.na(test)] <- 0
test <- test[,2:7]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/test[i,1]
}
test <- test[,-1]
rownames(test) <- my_level
write.csv(x = test,file = paste0("./output/1011/CD19_CD21_CD81_CD45freq.csv"))
test <- as.data.frame(t(test))

heatmap(as.matrix(t(test)),scale = "row",Colv = NA, Rowv = NA)
heatmap.2(as.matrix(t(test)), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c( "white", "firebrick3"))(50))

test2 <- as.data.frame(scales::percent(as.matrix(test), 0.01))
My_percnet <- function(x){return(paste(round(100*x, 2), "%", sep=""))}
test2 <- as.character(as.data.frame(lapply(test,formattable::percent,2)))
test2 <- as.data.frame(lapply(test,scales::percent,2))


# creating a function to compute percentage
my_percent <- function(num, digits = 2, ...) {      
  percentage <-formatC(num * 100, format = "f", digits = digits, ...)
  
  # appending "%" symbol at the end of
  # calculate percentage value
  paste0(percentage, "%")
}
test2 <- as.data.frame(lapply(test,my_percent,2))
rownames(test2) <- c("CD19","CD21","CD81","CD19+CD21+CD81","CD45")
pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="row",number_color = "black")
library("pheatmap")
# single cell plot
select_CD19 <- subset(x = pbmc, slot="data",subset = CD19 > 0.0)
select_CR2 <- subset(x = pbmc, slot="data",subset = CR2 > 0.0)
select_CD81 <- subset(x = pbmc, slot="data",subset = CD81 > 0.0)
select_all <- subset(x = pbmc, slot="data",subset = CD19 > 0.0 & CR2 > 0.0 & CD81 > 0.0)
VlnPlot(object = pbmc, features = "CD8A")


as.data.frame(lapply(test,scales::percent,2))

# single cell plot
select_CD19 <- subset(x = pbmc, slot="data",subset = CD19 > 0.0)
select_CR2 <- subset(x = pbmc, slot="data",subset = CR2 > 0.0)
select_CD81 <- subset(x = pbmc, slot="data",subset = CD81 > 0.0)
select_all <- subset(x = pbmc, slot="data",subset = CD19 > 0.0 & CR2 > 0.0 & CD81 > 0.0)
VlnPlot(object = pbmc, features = "CD8A")


as.data.frame(lapply(test,scales::percent,2))

####################################
## CD19 21 81 
####################################
features <- c("CD19","CR2","CD81")
FeaturePlot(pbmc,features = features ,label=T,cols=c("lightgrey","#DE1F1F"),label.size = 4,pt.size = 0.2)
ggsave("./output/1028/FeaturePlot.png",width = 13,height = 13)
FeaturePlot(pbmc,features = features[3] ,label=T,cols=c("lightgrey","#DE1F1F"),label.size = 4,pt.size = 0.2)
ggsave("./output/1028/FeaturePlot3.png",width = 8,height = 8)

DotPlot(pbmc, features = features) + RotatedAxis()+scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
ggsave("./output/1028/DotPlot.png",width = 6,height = 8)
VlnPlot(object = pbmc, features = features)
ggsave("./output/1028/VlnPlot.png",width = 15,height = 10)
VlnPlot(object = pbmc, features = features,stack = T,flip=T)
ggsave("./output/1028/stackedVlnPlot.png",width = 12,height = 5)


# CD19 CR2 CD81 table
select_CD19 <- subset(x = pbmc, slot="data",subset = CD19 > 0.0)
select_CR2 <- subset(x = pbmc, slot="data",subset = CR2 > 0.0)
select_CD81 <- subset(x = pbmc, slot="data",subset = CD81 > 0.0)
select_all <- subset(x = pbmc, slot="data",subset = CD19 > 0.0 & CR2 > 0.0 & CD81 > 0.0)

#select_PTPRC <- subset(x = pbmc, slot="data",subset = PTPRC > 0.0)
#PTPRC_num <- data.frame(table(Idents(select_PTPRC)))
cell_num <- data.frame(x=table(pbmc$cell_type))
colnames(cell_num)[1] <- "Var1"

P_19 <- as.data.frame(table(Idents(select_CD19)))
P_21 <- as.data.frame(table(Idents(select_CR2)))
P_81 <- as.data.frame(table(Idents(select_CD81)))
P_all <- as.data.frame(table(Idents(select_all)))

cell <- as.data.frame(my_level)
colnames(cell) <- "Var1"
test <- left_join(cell,cell_num,by="Var1")
test <- left_join(test,P_19,by="Var1")
test <- left_join(test,P_21,by="Var1")
test <- left_join(test,P_81,by="Var1")
test <- left_join(test,P_all,by="Var1")
#test <- left_join(test,PTPRC_num,by="Var1")

colnames(test) <- c("cell","cell_num","CD19","CD21","CD81","CD19+CD21+CD81")
test[is.na(test)] <- 0
write.csv(x = test,file = paste0("./output/1028/CD19_CD21_CD81count.csv"))

test <- test[,2:ncol(test)]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/test[i,1]
}
test <- test[,-1]
rownames(test) <- my_level
write.csv(x = test,file = paste0("./output/1028/CD19_CD21_CD81freq.csv"))
test <- as.data.frame(t(test))

# creating a function to compute percentage
my_percent <- function(num, digits = 2, ...) {      
  percentage <-formatC(num * 100, format = "f", digits = digits, ...)
  
  # appending "%" symbol at the end of
  # calculate percentage value
  paste0(percentage, "%")
}
test2 <- as.data.frame(lapply(test,my_percent,2))
rownames(test2) <- c("CD19","CD21","CD81","CD19+CD21+CD81")
write.csv(x = test2,file = paste0("./output/1028/CD19_CD21_CD81freq_percent.csv"))

pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="row",number_color = "black",fontsize = 18,gaps_row = c(3))
# png 1500 500 pdf 15 5

pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="column",number_color = "black",fontsize = 18,gaps_row = c(3))


####################################
## TCR 
####################################
# load TCR gene symbol name
TCR <- read.table("./data/TCR.txt")
features <- TCR$V1
# Umap plot
i=1
while (i <= length(features)){
  a = i
  b = i+1
  features1 <- features[a:b]
  print(features1)
  FeaturePlot(pbmc,features = features1 ,label=pbmc$cell_type,cols=c("lightgrey","#DE1F1F"),label.size = 5,pt.size = 0.2)
  ggsave(paste0("./output/1017/","FeaturePlot",{i},".pdf"),width = 20,height = 10)
  ggsave(paste0("./output/1017/","FeaturePlot",{i},".png"),width = 20,height = 10)
  i=i+2
}

# dotplot
{
  DotPlot(pbmc, features = features) + RotatedAxis()+scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
  ggsave("./output/1017/DotPlot.png",width = 8,height = 10)
  ggsave("./output/1017/DotPlot.pdf",width = 8,height = 10)
}

# violin plot
VlnPlot(object = pbmc, features = features)
VlnPlot(object = pbmc, features = features,stack = T,flip=T)
ggsave("./output/1017/VlnPlot.png",width = 15,height = 12)

## proportion table
# TCR table
features
#"CD247" "CD3D"  "CD3E"  "CD3G"  "TRAC"  "TRBC1" "TRBC2" "TRDC"  "TRGC1" "TRGC2"
{
  select_1 <- subset(x = pbmc, slot="data",subset = CD247 > 0.0)
  select_2 <- subset(x = pbmc, slot="data",subset = CD3D > 0.0)
  select_3 <- subset(x = pbmc, slot="data",subset = CD3E > 0.0)
  select_4 <- subset(x = pbmc, slot="data",subset = CD3G > 0.0)
  select_5 <- subset(x = pbmc, slot="data",subset = TRAC > 0.0)
  select_6 <- subset(x = pbmc, slot="data",subset = TRBC1 > 0.0)
  select_7 <- subset(x = pbmc, slot="data",subset = TRBC2 > 0.0)
  select_8 <- subset(x = pbmc, slot="data",subset = TRDC > 0.0)
  select_9 <- subset(x = pbmc, slot="data",subset = TRGC1 > 0.0)
  select_10 <- subset(x = pbmc, slot="data",subset = TRGC2 > 0.0)
  select_CD3 <- subset(x = pbmc, slot="data",subset = CD247 > 0.0 & CD3D>0.0 & CD3E>0.0 & CD3G>0.0)
  select_TCR1 <- subset(x = pbmc, slot="data",subset = TRAC > 0.0 & (TRBC1 > 0.0 |TRBC2 >0.0))
  select_TCR2 <- subset(x = pbmc, slot="data",subset = TRDC > 0.0 & (TRGC1 > 0.0 |TRGC2 >0.0))
  select_CD3_TCR1 <- subset(x = pbmc, slot="data",subset = (CD247 > 0.0 & CD3D>0.0 & CD3E>0.0 & CD3G>0.0) & (TRAC > 0.0 & (TRBC1 > 0.0 |TRBC2 >0.0)))
  select_CD3_TCR2 <- subset(x = pbmc, slot="data",subset = (CD247 > 0.0 & CD3D>0.0 & CD3E>0.0 & CD3G>0.0) & (TRDC > 0.0 & (TRGC1 > 0.0 |TRGC2 >0.0)))
}
{
  P_1 <- as.data.frame(table(Idents(select_1)))
  P_2 <- as.data.frame(table(Idents(select_2)))
  P_3 <- as.data.frame(table(Idents(select_3)))
  P_4 <- as.data.frame(table(Idents(select_4)))
  P_5 <- as.data.frame(table(Idents(select_5)))
  P_6<- as.data.frame(table(Idents(select_6)))
  P_7 <- as.data.frame(table(Idents(select_7)))
  P_8 <- as.data.frame(table(Idents(select_8)))
  P_9 <- as.data.frame(table(Idents(select_9)))
  P_10 <- as.data.frame(table(Idents(select_10)))
  P_CD3 <- as.data.frame(table(Idents(select_CD3)))
  P_TCR1 <- as.data.frame(table(Idents(select_TCR1)))
  P_TCR2 <- as.data.frame(table(Idents(select_TCR2)))
  P_CD3_TCR1 <- as.data.frame(table(Idents(select_CD3_TCR1)))
  P_CD3_TCR2 <- as.data.frame(table(Idents(select_CD3_TCR2)))
}

cell <- as.data.frame(my_level)
colnames(cell) <- "Var1"

cell_num <- data.frame(x=table(pbmc$cell_type))
colnames(cell_num)[1] <- "Var1"

test <- join_all(list(cell,cell_num,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_CD3,P_TCR1,P_TCR2,P_CD3_TCR1,P_CD3_TCR2), by='Var1', type='left')
test[is.na(test)] <- 0
colnames(test) <- c("cell","cell_num",features,"CD3 complex","TCRαβ","TCRγδ","TCRαβ complex","TCRγδ complex")

write.csv(x = test,file = paste0("./output/1017/TCR_count.csv"))

test <- test[,2:length(colnames(test))]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/test[i,1]
}
test <- test[,-1]
rownames(test) <- my_level
write.csv(x = test,file = paste0("./output/1017/TCRfreq.csv"))

test <- as.data.frame(t(test))
heatmap(as.matrix(t(test)),scale = "row",Colv = NA, Rowv = NA)
heatmap.2(as.matrix(t(test)), scale = "none", col = bluered(100), 
          trace = "none", density.info = "none")
pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = TRUE,
         color = colorRampPalette(c( "white", "firebrick3"))(50))

test2 <- as.data.frame(scales::percent(as.matrix(test), 0.01))
My_percnet <- function(x){return(paste(round(100*x, 2), "%", sep=""))}
test2 <- as.character(as.data.frame(lapply(test,formattable::percent,2)))
test2 <- as.data.frame(lapply(test,scales::percent,2))


# creating a function to compute percentage
my_percent <- function(num, digits = 2, ...) {      
  percentage <-formatC(num * 100, format = "f", digits = digits, ...)
  
  # appending "%" symbol at the end of
  # calculate percentage value
  paste0(percentage, "%")
}
test2 <- as.data.frame(lapply(test,my_percent,2))
rownames(test2) <- c(features,"CD3 complex","TCRαβ","TCRγδ","TCRαβ complex","TCRγδ complex")
write.csv(x = test2,file = paste0("./output/1017/TCRfreq_percent.csv"))

pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="row",number_color = "black",fontsize = 18)
# png 1500 1200 pdf 15 12

pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="column",number_color = "black",fontsize = 18)



####################################
## BCR 
####################################
# load BCR gene symbol name
BCR <- read.table("./data/BCR.txt")
features <- BCR$V1
# Umap plot
i=1
while (i <= length(features)){
  a = i
  b = i+1
  features1 <- features[a:b]
  print(features1)
  FeaturePlot(pbmc,features = features1 ,label=pbmc$cell_type,cols=c("lightgrey","#DE1F1F"),label.size = 5,pt.size = 0.2)
  ggsave(paste0("./output/1021/","FeaturePlot",{i},".pdf"),width = 20,height = 10)
  ggsave(paste0("./output/1021/","FeaturePlot",{i},".png"),width = 20,height = 10)
  i=i+2
}

# dotplot
{
  DotPlot(pbmc, features = features) + RotatedAxis()+scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
  ggsave("./output/1021/DotPlot.png",width = 8,height = 10)
  ggsave("./output/1021/DotPlot.pdf",width = 8,height = 10)
}

# violin plot
#VlnPlot(object = pbmc, features = features)
#ggsave("./output/1021/VlnPlot.png",width = 15,height = 12)
features <- features[-13]
i=1
while (i <= length(features)){
  a = i
  b = i+2
  features1 <- features[a:b]
  print(features1)
  VlnPlot(object = pbmc, features = features1,ncol=3)
  ggsave(paste0("./output/1021/","VlnPlot",{i},".pdf"),width = 15,height = 5)
  ggsave(paste0("./output/1021/","VlnPlot",{i},".png"),width = 15,height = 5)
  i=i+3
}

VlnPlot(object = pbmc, features = features,stack = T,flip=T)
ggsave("./output/1021/Stacked_VlnPlot.png",width = 15,height = 12)

## proportion table
# BCR table
features <- BCR$V1
features <- features[-13]
features
# "CD79A" "CD79B" "IGHA1" "IGHA2" "IGHD"  "IGHE"  "IGHG1" "IGHG2" "IGHG3" "IGHG4"
# "IGHM"  "IGKC"  "IGLC1" "IGLC2" "IGLC3" "IGLC7"
{
  select_1 <- subset(x = pbmc, slot="data",subset = CD79A > 0.0)
  select_2 <- subset(x = pbmc, slot="data",subset = CD79B > 0.0)
  select_3 <- subset(x = pbmc, slot="data",subset = IGHA1 > 0.0)
  select_4 <- subset(x = pbmc, slot="data",subset = IGHA2 > 0.0)
  select_5 <- subset(x = pbmc, slot="data",subset = IGHD > 0.0)
  select_6 <- subset(x = pbmc, slot="data",subset = IGHE > 0.0)
  select_7 <- subset(x = pbmc, slot="data",subset = IGHG1 > 0.0)
  select_8 <- subset(x = pbmc, slot="data",subset = IGHG2 > 0.0)
  select_9 <- subset(x = pbmc, slot="data",subset = IGHG3 > 0.0)
  select_10 <- subset(x = pbmc, slot="data",subset = IGHG4 > 0.0)
  select_11 <- subset(x = pbmc, slot="data",subset =  IGHM > 0.0)
  select_12 <- subset(x = pbmc, slot="data",subset =  IGKC> 0.0)
  #select_13 <- subset(x = pbmc, slot="data",subset =  IGLC1 > 0.0)
  select_14 <- subset(x = pbmc, slot="data",subset = IGLC2 > 0.0)
  select_15 <- subset(x = pbmc, slot="data",subset = IGLC3 > 0.0)
  select_16 <- subset(x = pbmc, slot="data",subset = IGLC7 > 0.0)
  select_79AB <- subset(x = pbmc, slot="data",subset = CD79A > 0.0 & CD79B >0.0)
  select_LC <- subset(x = pbmc, slot="data",subset =  (IGKC > 0.0 | IGLC1>0 |IGLC2>0.0|IGLC3>0.0|IGLC7>0.0))
  select_HC <- subset(x = pbmc, slot="data",subset = (IGHM > 0.0 | IGHG1>0.0|IGHG2>0.0|IGHG3>0.0|IGHG4>0.0|IGHA1>0.0|IGHA2>0|IGHE>0|IGHD>0))
  select_HL <- subset(x = pbmc, slot="data",subset = (IGKC > 0.0 | IGLC1>0 |IGLC2>0.0|IGLC3>0.0|IGLC7>0.0)&(IGHM > 0.0 | IGHG1>0.0|IGHG2>0.0|IGHG3>0.0|IGHG4>0.0|IGHA1>0.0|IGHA2>0|IGHE>0|IGHD>0))
  select_BCR <- subset(x = pbmc, slot="data",subset = (IGKC > 0.0 | IGLC1>0 |IGLC2>0.0|IGLC3>0.0|IGLC7>0.0)&(IGHM > 0.0 | IGHG1>0.0|IGHG2>0.0|IGHG3>0.0|IGHG4>0.0|IGHA1>0.0|IGHA2>0|IGHE>0|IGHD>0)&(CD79A > 0.0 & CD79B >0.0))
}
{
  P_1 <- as.data.frame(table(Idents(select_1)))
  P_2 <- as.data.frame(table(Idents(select_2)))
  P_3 <- as.data.frame(table(Idents(select_3)))
  P_4 <- as.data.frame(table(Idents(select_4)))
  P_5 <- as.data.frame(table(Idents(select_5)))
  P_6<- as.data.frame(table(Idents(select_6)))
  P_7 <- as.data.frame(table(Idents(select_7)))
  P_8 <- as.data.frame(table(Idents(select_8)))
  P_9 <- as.data.frame(table(Idents(select_9)))
  P_10 <- as.data.frame(table(Idents(select_10)))
  P_11 <- as.data.frame(table(Idents(select_11)))
  P_12 <- as.data.frame(table(Idents(select_12)))
  #P_13 <- as.data.frame(table(Idents(select_13))) IGLC1
  P_14 <- as.data.frame(table(Idents(select_14)))
  P_15 <- as.data.frame(table(Idents(select_15)))
  P_16 <- as.data.frame(table(Idents(select_16)))
  P_79AB <- as.data.frame(table(Idents(select_79AB)))
  P_LC <- as.data.frame(table(Idents(select_LC)))
  P_HC <- as.data.frame(table(Idents(select_HC)))
  P_HL <- as.data.frame(table(Idents(select_HL)))
  P_BCR <- as.data.frame(table(Idents(select_BCR)))
}

cell <- as.data.frame(my_level)
colnames(cell) <- "Var1"

cell_num <- data.frame(x=table(pbmc$cell_type))
colnames(cell_num)[1] <- "Var1"

#test <- join_all(list(cell,cell_num,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_13,P_14,P_15,P_16,P_79AB,P_HL,P_LC,P_HL,P_BCR), by='Var1', type='left')
test <- join_all(list(cell,cell_num,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_14,P_15,P_16,P_79AB,P_HL,P_LC,P_HL,P_BCR), by='Var1', type='left')
test[is.na(test)] <- 0
colnames(test) <- c("cell","cell_num",features,"CD79A+CD79B","Heavy Chain","Lignt Chain","Heavy Chain+Lignt Chain","BCR complex")
write.csv(x = test,file = paste0("./output/1021/BCR_count.csv"))

test <- test[,2:length(colnames(test))]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/test[i,1]
}
test <- test[,-1]
rownames(test) <- my_level
write.csv(x = test,file = paste0("./output/1021/BCRfreq.csv"))

#test$IGLC1 <- c(rep(0,nrow(test)))
#colnames(test)
#test <- test[,c(1:12,21,13:20)]
test <- as.data.frame(t(test))
# creating a function to compute percentage
my_percent <- function(num, digits = 2, ...) {      
  percentage <-formatC(num * 100, format = "f", digits = digits, ...)
  
  # appending "%" symbol at the end of
  # calculate percentage value
  paste0(percentage, "%")
}
test2 <- as.data.frame(lapply(test,my_percent,2))
rownames(test2) <- c(features,"CD79A+CD79B","Heavy Chain","Lignt Chain","Heavy Chain+Lignt Chain","BCR complex")
write.csv(x = test2,file = paste0("./output/1021/BCRfreq_percent.csv"))

save_pheatmap_pdf <- function(x, filename, width=15, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=1500, height=1200) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

p1 <- pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="row",number_color = "black",fontsize = 16,gaps_row = c(15))
save_pheatmap_pdf(p1,"./output/1021/heatmap_row.pdf",width = 15,height = 12)

p1 <- pheatmap(as.matrix(test),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2,
               scale="row",number_color = "black",fontsize = 18,gaps_row = c(15))
save_pheatmap_png(p1,"./output/1021/heatmap_row.png",width = 1200,height = 1000)


pheatmap(as.matrix(test),
         cluster_row = FALSE,cluster_cols = FALSE,
         display_numbers = test2,
         scale="column",number_color = "black",fontsize = 16,gaps_row = c(15))


p2 <- pheatmap(as.matrix(test),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2,
               scale="column",number_color = "black",fontsize = 16,gaps_row = c(15))
save_pheatmap_pdf(p2,"./output/1021/heatmap_col.pdf",width = 15,height = 12)

p2 <- pheatmap(as.matrix(test),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2,
               scale="column",number_color = "black",fontsize = 18,gaps_row = c(15))
save_pheatmap_png(p2,"./output/1021/heatmap_col.png",width = 1200,height = 1000)

#############################################
#CD4 CD8
select_CD4 <- subset(x = pbmc, slot="data",subset = CD4 > 0.0)
select_CD8A <- subset(x = pbmc, slot="data",subset = CD8A > 0.0)

P_4 <- as.data.frame(table(Idents(select_CD4)))
P_8 <- as.data.frame(table(Idents(select_CD8A)))
test <- left_join(cell_num,P_4,by="Var1")
test <- left_join(test,P_8,by="Var1")
colnames(test) <- c("cell","cell_num","CD4","CD8")
write.csv(x = test,file = paste0("./output/1006/CD4_CD8.csv"))
test[is.na(test)] <- 0
test <- test[,3:4]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/cell_num[i,2]
}
write.csv(x = test,file = paste0("./output/1006/CD14_CD8_freq.csv"))



Human_PBMC<-readRDS("./data/blish_covid.seu.rds")
DefaultAssay(Human_PBMC)<-"RNA"
Idents(Human_PBMC) <- Human_PBMC$Status
table(Idents(Human_PBMC))
pbmc<-subset(x = Human_PBMC, idents = "Healthy") #16627 cell 6healthy
Idents(pbmc) <- pbmc$cell.type.fine
features <- c("CD19","PSMB5","CR2","CD81")
FeaturePlot(pbmc,features = features ,label=pbmc$cell.type.fine,cols=c("lightgrey","#DE1F1F"),label.size = 4,pt.size = 0.2)
FeaturePlot(pbmc,features = features,cols=c("lightgrey","#DE1F1F"),label.size = 3.5,pt.size = 0.2,label = T)
DotPlot(pbmc,features = features,cols="Spectral",cluster.idents=TRUE) +theme(axis.text.x = element_text(hjust = 1,angle = 60)) 
FeaturePlot(pbmc,features = features, label = T, cols=c("lightgrey","#DE1F1F"),label.size = 4,pt.size = 0.2)


#######################
## MHC complex
#######################

# load MHC gene symbol name
MHC <- read.table("./data/MHC.txt")
features <- MHC$V1
features
features <- features[-c(20,21)]
# Umap plot
i=1
while (i <= length(features)){
  a = i
  b = i+1
  features1 <- features[a:b]
  print(features1)
  FeaturePlot(pbmc,features = features1 ,label=pbmc$cell_type,cols=c("lightgrey","#DE1F1F"),label.size = 5,pt.size = 0.2)
  ggsave(paste0("./output/1110/","FeaturePlot",{i},".pdf"),width = 20,height = 10)
  ggsave(paste0("./output/1110/","FeaturePlot",{i},".png"),width = 20,height = 10)
  i=i+2
}

# dotplot
{
  DotPlot(pbmc, features = features) + RotatedAxis()+scale_color_gradientn(colours = rev(brewer.pal(n = 5, name ="RdYlBu")))
  ggsave("./output/1110/DotPlot.png",width = 10,height = 10)
  ggsave("./output/1110/DotPlot.pdf",width = 10,height = 10)
}

# violin plot
#VlnPlot(object = pbmc, features = features)
#ggsave("./output/1021/VlnPlot.png",width = 15,height = 12)

i=1
while (i <= length(features)){
  a = i
  b = i+2
  features1 <- features[a:b]
  print(features1)
  VlnPlot(object = pbmc, features = features1,ncol=3)
  ggsave(paste0("./output/1110/","VlnPlot",{i},".pdf"),width = 15,height = 5)
  ggsave(paste0("./output/1110/","VlnPlot",{i},".png"),width = 15,height = 5)
  i=i+3
}

VlnPlot(object = pbmc, features = features,stack = T,flip=T)
ggsave("./output/1110/Stacked_VlnPlot.png",width = 15,height = 12)

## proportion table
# MHC table
features <- MHC$V1
features <- features[-c(20,21)]
features
# "HLA-A"    "HLA-B"    "HLA-C"    "HLA-E"    "HLA-F"    "HLA-G"    "B2M"     
# "HLA-DMA"  "HLA-DMB"  "HLA-DOA"  "HLA-DOB"  "HLA-DPA1" "HLA-DPB1" "HLA-DQA1"
# "HLA-DQA2" "HLA-DQB1" "HLA-DQB2" "HLA-DRA"  "HLA-DRB1" "HLA-DRB5"
{
  select_1 <- subset(x = pbmc, slot="data",subset = `HLA-A` > 0.0)
  select_2 <- subset(x = pbmc, slot="data",subset = `HLA-B` > 0.0)
  select_3 <- subset(x = pbmc, slot="data",subset = `HLA-C` > 0.0)
  select_4 <- subset(x = pbmc, slot="data",subset = `HLA-E` > 0.0)
  select_5 <- subset(x = pbmc, slot="data",subset = `HLA-F` > 0.0)
  select_6 <- subset(x = pbmc, slot="data",subset = `HLA-G` > 0.0)
  select_7 <- subset(x = pbmc, slot="data",subset = `HLA-DMA` > 0.0)
  select_8 <- subset(x = pbmc, slot="data",subset = `HLA-DMB` > 0.0)
  select_9 <- subset(x = pbmc, slot="data",subset = `HLA-DOA` > 0.0)
  select_10 <- subset(x = pbmc, slot="data",subset = `HLA-DOB` > 0.0)
  select_11 <- subset(x = pbmc, slot="data",subset =  `HLA-DPA1` > 0.0)
  select_12 <- subset(x = pbmc, slot="data",subset =  `HLA-DPB1` > 0.0)
  select_13 <- subset(x = pbmc, slot="data",subset =  `HLA-DQA1` > 0.0)
  select_14 <- subset(x = pbmc, slot="data",subset = `HLA-DQA2` > 0.0)
  select_15 <- subset(x = pbmc, slot="data",subset = `HLA-DQB1` > 0.0)
  select_16 <- subset(x = pbmc, slot="data",subset = `HLA-DQB2` > 0.0)
  select_17 <- subset(x = pbmc, slot="data",subset = `HLA-DRA` > 0.0)
  select_18 <- subset(x = pbmc, slot="data",subset = `HLA-DRB1` > 0.0)
  select_19 <- subset(x = pbmc, slot="data",subset = `HLA-DRB5` > 0.0)
  select_20 <- subset(x = pbmc, slot="data",subset = B2M > 0.0)
  #select_MHC1_a <- subset(x = pbmc, slot="data",subset = (`HLA-A` > 0.0 | `HLA-B` > 0.0 | `HLA-C` > 0.0 & `HLA-E` > 0.0 & `HLA-F` > 0.0 & `HLA-G` > 0.0))
  select_MHC1 <- subset(x = pbmc, slot="data",subset = ((`HLA-A` > 0.0 & B2M >0) | (`HLA-B` > 0.0 & B2M >0)|(`HLA-C` > 0.0 & B2M >0)|(`HLA-E` > 0.0 & B2M >0)|(`HLA-F` > 0.0 & B2M >0)|(`HLA-G` > 0.0 & B2M >0)))
  #select_classicMHC1 <- subset(x = pbmc, slot="data",subset = (`HLA-A` > 0.0 | `HLA-B` > 0.0 | `HLA-C` > 0.0))
  select_MHC2 <- subset(x = pbmc, slot="data",subset = (`HLA-DMA` > 0.0 & `HLA-DMB` > 0.0)|(`HLA-DOA` > 0.0 & `HLA-DOB` > 0.0)|(`HLA-DPA1` > 0.0& `HLA-DPB1`> 0.0)|(`HLA-DQA1` > 0.0& `HLA-DQB1`> 0.0)|(`HLA-DQA2` > 0.0& `HLA-DQB2`> 0.0) |(`HLA-DRA` > 0.0& `HLA-DRB1`> 0.0)|(`HLA-DRA` > 0.0& `HLA-DRB5`> 0.0))
  #select_classicMHC2 <- subset(x = pbmc, slot="data",subset = (`HLA-DPA1` > 0.0& `HLA-DPB1`> 0.0)|(`HLA-DQA1` > 0.0 & `HLA-DQB1`> 0.0)|(`HLA-DQA2` > 0.0& `HLA-DQB2`> 0.0) |(`HLA-DRA` > 0.0& `HLA-DRB1`> 0.0)|(`HLA-DRA` > 0.0& `HLA-DRB5`> 0.0))
}
{
  P_1 <- as.data.frame(table(Idents(select_1)))
  P_2 <- as.data.frame(table(Idents(select_2)))
  P_3 <- as.data.frame(table(Idents(select_3)))
  P_4 <- as.data.frame(table(Idents(select_4)))
  P_5 <- as.data.frame(table(Idents(select_5)))
  P_6<- as.data.frame(table(Idents(select_6)))
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
  P_17 <- as.data.frame(table(Idents(select_17)))
  P_18 <- as.data.frame(table(Idents(select_18)))
  P_19 <- as.data.frame(table(Idents(select_19)))
  P_20 <- as.data.frame(table(Idents(select_20)))
  P_MHC1 <- as.data.frame(table(Idents(select_MHC1)))
  P_MHC2 <- as.data.frame(table(Idents(select_MHC2)))
  #P_classicMHC1 <- as.data.frame(table(Idents(select_classicMHC1)))
  #P_classicMHC2 <- as.data.frame(table(Idents(select_classicMHC2)))
}

cell <- as.data.frame(my_level)
colnames(cell) <- "Var1"

cell_num <- data.frame(x=table(pbmc$cell_type))
colnames(cell_num)[1] <- "Var1"

test <- join_all(list(cell,cell_num,P_1,P_2,P_3,P_4,P_5,P_6,P_20,P_7,P_8,P_9,P_10,P_11,P_12,P_13,P_14,P_15,P_16,P_17,P_18,P_19,P_MHC1,P_MHC2), by='Var1', type='left')
#test <- join_all(list(cell,cell_num,P_1,P_2,P_3,P_4,P_5,P_6,P_7,P_8,P_9,P_10,P_11,P_12,P_14,P_15,P_16,P_79AB,P_HL,P_LC,P_HL,P_BCR), by='Var1', type='left')
test[is.na(test)] <- 0
colnames(test) <- c("cell","cell_num",features,"MHC I","MHC II")
write.csv(x = test,file = paste0("./output/1110/MHC_count.csv"))

test <- test[,2:length(colnames(test))]
for (i in 1:nrow(test)){
  test[i,] <- test[i,]/test[i,1]
}
test <- test[,-1]
rownames(test) <- my_level
write.csv(x = test,file = paste0("./output/1110/MHCfreq.csv"))

#test$IGLC1 <- c(rep(0,nrow(test)))
#colnames(test)
#test <- test[,c(1:12,21,13:20)]
test <- as.data.frame(t(test))
# creating a function to compute percentage
my_percent <- function(num, digits = 2, ...) {      
  percentage <-formatC(num * 100, format = "f", digits = digits, ...)
  
  # appending "%" symbol at the end of
  # calculate percentage value
  paste0(percentage, "%")
}
test2 <- as.data.frame(lapply(test,my_percent,2))
rownames(test2) <- c(features,"MHC I","MHC II")
write.csv(x = test2,file = paste0("./output/1110/MHCfreq_percent.csv"))

save_pheatmap_pdf <- function(x, filename, width=15, height=12) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=1500, height=1200) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  png(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

## single
p1 <- pheatmap(as.matrix(test[1:20,]),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2[1:20,],
               scale="row",number_color = "black",fontsize = 16)
save_pheatmap_pdf(p1,"./output/1110/single_heatmap_row.pdf",width = 15,height = 12)

p2 <- pheatmap(as.matrix(test[1:20,]),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2[1:20,],
               scale="column",number_color = "black",fontsize = 16)
save_pheatmap_pdf(p2,"./output/1110/single_heatmap_col.pdf",width = 15,height = 12)

## complex
p1 <- pheatmap(as.matrix(test[21:22,]),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2[21:22,],
               scale="row",number_color = "black",fontsize = 16)
save_pheatmap_pdf(p1,"./output/1110/complex_heatmap_row.pdf",width = 15,height = 3.3)

p2 <- pheatmap(as.matrix(test[21:22,]),
               cluster_row = FALSE,cluster_cols = FALSE,
               display_numbers = test2[21:22,],
               scale="column",number_color = "black",fontsize = 16)
save_pheatmap_pdf(p2,"./output/1110/complex_heatmap_col.pdf",width = 15,height = 3.3)







