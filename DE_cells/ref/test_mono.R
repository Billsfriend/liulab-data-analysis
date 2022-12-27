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
  library(plyr)
  library(pheatmap)
  options(stringsAsFactors = F)
}


#20211115
Human_PBMC<-readRDS("./data/Final_nCoV_0716_upload.RDS")
DefaultAssay(Human_PBMC)<-"RNA"
pbmc<-subset(x = Human_PBMC, idents = "Ctrl")

table(pbmc$cell_type)
table(Idents(pbmc))
#all.genes <- rownames(pbmc) 
all.genes <- rownames(pbmc@assays$RNA@data)
my_level <- c("Naive B cells","Memory B cells","Plasma","Cycling Plasma","Naive T cells","Activated CD4 T cells","Cytotoxic CD8 T cells","Cycling T cells","MAIT","NKs","XCL+ NKs","Monocytes","Megakaryocytes","DCs","Stem cells")
Idents(pbmc) <- factor(pbmc$cell_type,levels = my_level)

monocytes <- subset(x = pbmc, slot="data",idents = "Monocytes")
mono_CD79B <- subset(x = pbmc, slot="data",subset = CD79B > 0, idents = "Monocytes")
mono_CD79B_2 <- subset(x = pbmc, slot="data",subset = CD79B == 0, idents = "Monocytes")
Idents(object = monocytes, cells = colnames(mono_CD79B)) <- "CD79B+Monocytes"
Idents(object = monocytes, cells = colnames(mono_CD79B_2)) <- "CD79B-Monocytes"

#mono_CD79B.marker <- FindMarkers(monocytes, ident.1 = "mono_CD79B", min.pct = 0.1)
#mono_CD79B.marker <- FindMarkers(monocytes, ident.1 = "mono_CD79B", min.pct = 0.25)
#mono_CD79B.marker["FCGR3A",]
#table(Idents(monocytes))
mono_CD79B.marker <- FindMarkers(monocytes, ident.1 = "CD79B+Monocytes", ident.2= "CD79B-Monocytes",min.pct = 0.1,test.use = "MAST")
which(rownames(mono_CD79B.marker)=="FCGR3A")
#write.table(mono_CD79B.marker,"./output/1115/BGI-mono-gene.txt",sep = "\t",quote = F)
#visualization
genelist <-  rownames(mono_CD79B.marker)
#VlnPlot(object = monocytes, features =genelist[1:5])

BGI.mono <- mono_CD79B.marker
#count cell number
mono_CD79B_CD16 <- subset(x = monocytes, slot="data",subset = CD79B> 0 & FCGR3A> 0) #39
mono_CD14 <- subset(x = monocytes, slot="data",subset = CD14 > 0) #6
mono_FCGR3A <- subset(x = monocytes, slot="data",subset = FCGR3A > 0) #104
monocytes # 132
mono_CD79B <- subset(x = monocytes, slot="data",subset = CD79B>0) #40

monocytes@assays$RNA@scale.data <- scale(monocytes@assays$RNA@data, scale = TRUE)
?DoHeatmap
DoHeatmap(monocytes,features = genelist[1:30],group.by = "ident",angle=0)


# -------------------------------------------
# Seq-Well Data
# -------------------------------------------

Human_PBMC<-readRDS("./data/blish_covid.seu.rds")
DefaultAssay(Human_PBMC)<-"RNA"
Idents(Human_PBMC) <- Human_PBMC$Status
table(Idents(Human_PBMC))
pbmc<-subset(x = Human_PBMC, idents = "Healthy") #16627 cell 6healthy
Idents(pbmc) <- pbmc$cell.type.fine
x <- as.data.frame(table(Idents(pbmc)))
#Idents(Human_PBMC) <- Human_PBMC$cell.type.fine
y <- as.data.frame(table(Idents(Human_PBMC)))

c("RBC","CD4n T","CD14 Monocyte","B","CD8m T","NK","CD16 Monocyte","SC & Eosinophil","Platelet","Neutrophil","IgA PB","pDC","DC" ,"CD4m T","CD8eff T",
  "gd T","CD4 T","IgG PB","Activated Granulocyte","Class-switched B")
my_level <- c("B","Class-switched B","IgA PB","IgG PB","CD4 T","gd T","CD4n T","CD4m T","CD8eff T","CD8m T","NK","CD14 Monocyte","CD16 Monocyte","pDC","DC", "Neutrophil","SC & Eosinophil","Activated Granulocyte","Platelet","RBC")
#my_level <- c("Naive B cells","Memory B cells","Plasma","Cycling Plasma","Naive T cells","Activated CD4 T cells","Cytotoxic CD8 T cells","Cycling T cells","MAIT","NKs","XCL+ NKs","Monocytes","Megakaryocytes","DCs","Stem cells")
Idents(pbmc) <- factor(pbmc$cell.type.fine,levels = my_level)
table(Idents(pbmc))


monocytes <- subset(x = pbmc, slot="data",idents = c("CD14 Monocyte","CD16 Monocyte"))
mono16_CD79B <- subset(x = pbmc, slot="data",subset = CD79B > 0, idents = "CD16 Monocyte")
mono14_CD79B <- subset(x = pbmc, slot="data",subset = CD79B > 0, idents = "CD14 Monocyte")
mono16_CD79B_2 <- subset(x = pbmc, slot="data",subset = CD79B == 0, idents = "CD16 Monocyte")
mono14_CD79B_2 <- subset(x = pbmc, slot="data",subset = CD79B == 0, idents = "CD14 Monocyte")
Idents(object = monocytes, cells = c(colnames(mono16_CD79B),colnames(mono14_CD79B))) <- "CD79B+Monocytes"
Idents(object = monocytes, cells = c(colnames(mono16_CD79B_2),colnames(mono14_CD79B_2))) <- "CD79B-Monocytes"
mono_CD79B.marker <- FindMarkers(monocytes, ident.1 = "CD79B+Monocytes", ident.2= "CD79B-Monocytes",min.pct = 0.1,test.use = "MAST")
which(rownames(mono_CD79B.marker)=="FCGR3A")
new.mono <- mono_CD79B.marker
#write.table(mono_CD79B.marker,"./output/1115/newdata-mono-gene.txt",sep = "\t",quote = F)
#visualization
genelist2 <-  rownames(mono_CD79B.marker)
VlnPlot(object = monocytes, features =genelist[1:5])

mono_CD79B.marker["FCGR3A",]

#火山图
library(ggpubr)
library(ggthemes)
library(ggrepel)

BGI <- read.table("./output/1115/BGI-mono-gene.txt",sep = "\t",row.names = 1)
new <- read.table("./output/1115/newdata-mono-gene.txt",sep = "\t",row.names = 1)

BGI <- BGI[BGI$p_val<0.05,] #437
new <- new[new$p_val<0.05,] #2006
BGI_pos <- BGI[BGI$avg_log2FC>0,] 
BGI_neg <- BGI[BGI$avg_log2FC<0,] 
new_pos <- new[new$avg_log2FC>0,]
new_neg <- new[new$avg_log2FC<0,]

# 对差异p(adj.P.Val一列进行log10转换
BGI$logP <- -log10(BGI$p_val)
new$logP <- -log10(new$p_val)

df3 <- new
df3 <- df3[-1,]
df3$threshold <- as.factor(ifelse(df3$p_val < 0.05 & abs(df3$avg_log2FC) >=1,ifelse(df3$avg_log2FC > 1 ,'Up','Down'),'Not'))
p <- ggplot(data=df3, aes(x=avg_log2FC, y =-log10(p_val), color=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(size=1) +
  xlim(c(-10, 10)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2FoldChange",y="-log10 (p-value)",title="new data", face="bold")
df3$symbol <- rownames(df3)
df3$label <- ifelse(df3$p_val < 0.05 & abs(df3$avg_log2FC) >=5,df3$symbol,"")

#p + geom_text_repel(aes(avg_log2FC, df3$avg_log2FC, label=label))

p+geom_text_repel(data = df3, aes(x = avg_log2FC, 
                                   y = -log10(p_val), 
                                   label = label),
                  size = 3,box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.8, "lines"), 
                  segment.color = "black", 
                  show.legend = FALSE,max.overlaps = 300)

p
#韦恩图 交集gene list
"FCGR3A" %in% genelist_pos
genelist <- intersect(rownames(BGI),rownames(new))
genelist_pos <- intersect(rownames(BGI_pos),rownames(new_pos)) #两个data 同时上调的
genelist_neg <- intersect(rownames(BGI_neg),rownames(new_neg)) #两个data 同时下调的

BGI <- BGI[order(BGI$avg_log2FC,decreasing = T),]
new <- new[order(new$avg_log2FC,decreasing = T),]

#write.table(BGI,"./output/1115/BGI-mono-gene-sort.txt",sep = "\t",quote = F)
#write.table(new,"./output/1115/New-mono-gene-sort.txt",sep = "\t",quote = F)

genelist_pos <- intersect(rownames(BGI_pos),rownames(new_pos)) #两个data 同时上调的
genelist_neg <- intersect(rownames(BGI_neg),rownames(new_neg)) #两个data 同时下调的
write.table(genelist_pos,"./output/1115/genelist_pos.txt")
write.table(genelist_neg,"./output/1115/genelist_neg.txt")
#提供gene list

#尝试共表达模块和heatmap
which(rownames(BGI)=="FCGR3A")
df3["FCGR3A",]

colnames(BGI_pos) <- paste0("BGI_",colnames(BGI_pos))
BGI_pos$name <- rownames(BGI_pos)
colnames(new_pos) <- paste0("Newdata_",colnames(new_pos))
new_pos$name <- rownames(new_pos)
pos.gene.info <- inner_join(BGI_pos,new_pos,by="name")
pos.gene.info <- pos.gene.info[,c(6,1:5,7:11)]

colnames(BGI_neg) <- paste0("BGI_",colnames(BGI_neg))
BGI_neg$name <- rownames(BGI_neg)
colnames(new_neg) <- paste0("Newdata_",colnames(new_neg))
new_neg$name <- rownames(new_neg)
neg.gene.info <- inner_join(BGI_neg,new_neg,by="name")
neg.gene.info <- neg.gene.info[,c(6,1:5,7:11)]

pos.genelist <- pos.gene.info[,1]
neg.genelist <- neg.gene.info[,1]

#write.table(pos.genelist,"./output/1121/pos_genelist.txt",quote = F, row.names = F,col.names = F)
#write.table(neg.genelist,"./output/1121/neg_genelist.txt",quote = F, row.names = F,col.names = F)

pos_pro <- read.table("output/1121/pos_gene_protein.txt",sep = "\t")
colnames(pos_pro) <- pos_pro[1,]
pos_pro <- pos_pro[-1,]
pos.gene.info$Uniprot_Protein_name <- pos_pro[,2]
pos.gene.info <- pos.gene.info[,c(1,12,2:11)]
write.table(pos.gene.info,"./output/1121/pos.gene.info.txt",quote = F,sep = "\t",row.names = F)
write.csv(pos.gene.info,"./output/1121/pos.gene.info.csv",quote = F,row.names = F)


neg_pro <- read.table("output/1121/neg_gene_protein.txt",sep = "\t")
colnames(neg_pro) <- neg_pro[1,]
neg_pro <- neg_pro[-1,]
neg.gene.info$Uniprot_Protein_name <- neg_pro[,2]
neg.gene.info <- neg.gene.info[,c(1,12,2:11)]
write.table(pos.gene.info,"./output/1121/neg.gene.info.txt",quote = F,sep = "\t",row.names = F)
write.csv(neg.gene.info,"./output/1121/neg.gene.info.csv",quote = F,row.names = F)


## -------------------------------------------
## funtional_analysis 
genelist <- neg.gene.info[,1]
geneid <-  bitr(genelist,fromType="SYMBOL",toType = c("ENTREZID"),OrgDb = "org.Hs.eg.db")

#KEGG
KEGG <- enrichKEGG(
  gene = geneid$ENTREZID, 
  keyType = "kegg", 
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  qvalueCutoff = 1
)
dotplot(KEGG,showCategory=12, font.size=15, title="KEGG Pathway Enrichment") 

GO<-enrichGO(geneid$ENTREZID,#GO富集分析
             OrgDb = "org.Hs.eg.db",
             keyType = "ENTREZID",#设定读取的gene ID类型
             ont = "ALL",#(ont为ALL因此包括 Biological Process,Cellular Component,Mollecular Function三部分）
             pvalueCutoff = 1,#设定p值阈值#设定q值阈值
             qvalueCutoff = 1,
             pAdjustMethod="BH",
             readable = T)
dotplot(GO, split="ONTOLOGY")+facet_grid(ONTOLOGY~., scale="free")
