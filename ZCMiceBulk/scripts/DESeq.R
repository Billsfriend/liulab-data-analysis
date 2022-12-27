# ref: https://www.jianshu.com/p/bdb871a08506
# BiocManager::install('DESeq2')
library(DESeq2)
library(pheatmap)
library(BiocParallel)
library(RColorBrewer)
library(sva)

# read count data
counts <- read.csv("data/zccounts.csv", row.names = 1)

# transform to numeric matrix
# Count <- data.matrix(counts)

# define group factor
zccoldata <- read.csv('DiffAnalysis/sampleinfo.txt', sep = '\t', row.names = 1)

# generate DESeqDataSet object
zcdds <- DESeqDataSetFromMatrix(counts, DataFrame(zccoldata), ~ condition)

# filter low expression data
zcdds <- zcdds[rowSums(counts(dds)) > 10, ] 

# relevel control group as reference
dds$condition <- relevel(dds$condition, ref = "control")

# regularized-logarithm transformation. transform counts data to log2
rld <- rlog(object=dds,blind=F) 
# or use variance stabilizing transformation instead of rlog:
vsd <- vst(object=dds,blind=T)

# estimate size factor of data. prepare for normalization
dds <- estimateSizeFactors(dds) 

# normalized counts = read counts/sizeFactor
normalized_counts <- counts(dds,normalized=T) 

# remove hidden batch effect
# dat <- counts(dds, normalized=TRUE)
# idx <- rowMeans(dat) > 1
# dat <- dat[idx,]
# mod <- model.matrix(~ dat, colData(dds))
# mod0 <- model.matrix(~ 1, colData(dds))
# svseq <- svaseq(dat, mod, mod0, n.sv=2)
#### encounter error in the line above
# svseq$sv
# par(mfrow=c(2,1),mar=c(3,5,3,1))
# stripchart(svseq$sv[,1] ~ dds$cell,vertical=TRUE,main="SV1")
# abline(h=0)
# stripchart(svseq$sv[,2] ~ dds$cell,vertical=TRUE,main="SV2")
# abline(h=0)
# ddssva <- dds
# ddssva$SV1 <- svseq$sv[,1]
# ddssva$SV2 <- svseq$sv[,2]
# design(ddssva) <- ~ SV1 + SV2 + dex
# ddssva <- DESeq(ddssva)

dds <- DESeq(dds)
# equal to following stepwise analysis
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds) 
# dds <- nbinomWaldTest(dds)
 
# obtain results of analysis
res <- results(dds)

### Output DE result

# order by adjusted p value, from low to high
resOrdered <- res[order(res$padj), ]
summary(resOrdered)

# transform DESeqResults object to dataframe for output
resdata <- merge(as.data.frame(resOrdered), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

names(resdata)[1] <- "Gene"

# select only genes of p.adj < 0.05 and log2FC > 1
resdata_filter <- subset(resdata,resdata$padj<0.05 & abs(resdata$log2FoldChange) >1)

# write DEG data to disk
write.table(resdata_filter, file="DiffAnalysis/diffexpr_Sample1_Sample2_0.05.txt",sep="\t",row.names=F)

write.table(resdata,file="DiffAnalysis/diffexpr_Sample1_Sample2_results.all.txt",quote=F,sep="\t",row.names=F)

# shrink logFC to correct its value
resultsNames(dds) 

# apeglm package is needed
resLFC <- lfcShrink(dds, coef="condition_knockout_vs_control", type="apeglm")

resShrink <- merge(as.data.frame(resLFC), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)

write.table(resShrink,file="DiffAnalysis/diffexpr_Sample1_Sample2_results.LFCshrink.txt",quote=F,sep="\t",row.names=F)

# use likelihood ratio test to perform DESeq
# ddsLRT <- DESeq(dds, test="LRT", reduced=~1)
# resLRT <- results(ddsLRT)


register(MulticoreParam(4))
#先预定4个核，等需要的时候直接使用parallel=TRUE来调用。

# visualize with MA plot

colData <- read.table("DiffAnalysis/sampleinfo.txt", header = T)

Pairindex <- read.table("DiffAnalysis/compinfo.txt", sep = "\t", header = T)

for(i in 1:nrow(Pairindex)){
  
  j<-as.vector(unlist(Pairindex[i,][1]))
  
  k<-as.vector(unlist(Pairindex[i,][2]))
  
  print(paste(j,k,sep=" vs "))
  
  sample_1 = j
  
  sample_2 = k
  
}

png(paste0("figures/",sample_1,"_vs_",sample_2,"_MAplot.png"))

plotMA(res, main=paste0(sample_1,"_vs_",sample_2,("_MAplot")), ylim=c(-10,10))

dev.off()

# draw heatmap of corelation between samples by pheatmap
ddCor <- cor(normalized_counts, method = "pearson")

pheatmap(file = "figures/control_VS_knockout_cor-zc.png", ddCor, clustering_method = "average", display_numbers = T)

# use rlog to draw heatmap
rld
sampleDist <- dist(t(assay(rld)))  #样品距离,欧氏距离,t转置 
#为确保所有基因大致相同的contribution用rlog-transformed data 
#画某些基因在样本间的heatmap也可以用rlog数据 
#用PoiClaClu包计算泊松距离(Poisson Distance)，必须是原始表达矩阵 
#poisd <- PoissonDistance(t(counts(dds))) 
sampleDistMatrix <- as.matrix(sampleDist)  #样品间距离的矩阵

head(sampleDistMatrix)  #样品间距离的数据框
colors <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDist,
         clustering_distance_cols=sampleDist,
         color = colors)

library(ggplot2)
#把样本之间的距离转化成二维坐标,在降维过程中保证样品点之间的距离不变
#MDS/PCoA基于最小线性距离(欧氏距离)的聚类,与PCA的最大线性相关是一样的
#适合在没有表达矩阵值,但只有一个距离矩阵的情况下使用
mdsdata <- data.frame(cmdscale(sampleDistMatrix))
#cmdscale(classical multidimensional scaling)
#返回MDS的坐标
mds <- cbind(mdsdata,as.data.frame(colData(vsd)))
mds  #按列合并
ggplot(data=mds,aes(X1,X2,color=Count_condition)) +
  geom_point(size=3)

# PCA
pcadata <- plotPCA(vsd, intgroup = 'condition', returnData=TRUE)
percentVar <- round(100*attr(pcadata,"percentVar"),1)
ggplot(pcadata, aes(PC1, PC2, color=condition)) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_bw()

topGene <- rownames(res)[which.min(res$padj)]  #padj最小的一个基因
plotCounts(dds, gene=topGene, intgroup=c("ko"))  #画出这个基因的标准化后的表达量
#以散点图的形式画出这个基因在各样本中的表达量
data <- plotCounts(dds, gene=topGene, intgroup=c("dex","cell"), returnData=TRUE)
ggplot(data, aes(x=dex, y=count, color=cell)) +  
  scale_y_log10() +
  geom_point(position=position_jitter(width=.1,height=0), size=3)  
ggplot(data, aes(x=dex, y=count, fill=dex)) +
  scale_y_log10() +
  geom_dotplot(binaxis="y", stackdir="center")
ggplot(data, aes(x=dex, y=count, color=cell, group=cell)) +
  scale_y_log10() + geom_point(size=3) + geom_line()

plotMA(res, ylim=c(-5,5))
#"M" for minus（减）, because a log ratio is equal to log minus log, and "A" for average（均值）
#M对应差异对比组之间基因表达变化log2 fold changes (Y轴)
#A对应差异对比组基因表达量均值the mean of normalized counts (X轴)
plotMA(res, alpha = 0.1, main = "", xlab = "mean of normalized counts", ylim=c(-5,5))
#alpha为padj显著性水平阈值，默认alpha=0.1
#Each gene is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
plotMA(resLFC, ylim=c(-5,5))  #提高了log2 fold change阈值
topGene <- rownames(resLFC)[which.min(resLFC$padj)]
with(resLFC[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})  #标记出一个特定的基因

# volcano plot
# g <- read.table("diffexpr_Sample1_Sample2_results.all.txt", header = T, row.names = 1)
g <- resdata

g <- na.omit(g)

plot(g$log2FoldChange, g$padj)

g <- transform(g, padj = -1*log10(g$padj))

down <- g[g$log2FoldChange <= -1 & g$padj > 3,]

up <- g[g$log2FoldChange >=1 & g$padj > 3,]

no <- g[g$log2FoldChange > -1 & g$log2FoldChange < 1 | g$padj < 3,]


png(filename = 'figures/volc0.001.png')
plot(no$log2FoldChange, no$padj, xlim = c(-10,10), ylim = c(0,100), col = "blue", pch = 16, cex = 0.8, main = "Gene Expression", xlab = "log2FoldChange", ylab = "-log10(Qvalue)")
points(up$log2FoldChange, up$padj, col = "red", pch = 16, cex = 0.8)
points(down$log2FoldChange, down$padj, col = "green", pch = 16, cex = 0.8)
dev.off()

