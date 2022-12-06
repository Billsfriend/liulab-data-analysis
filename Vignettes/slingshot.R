# BiocManager::install('slingshot')

library(slingshot)

browseVignettes('slingshot')

# generate synthetic count data representing a single lineage
means <- rbind(
  # non-DE genes
  matrix(rep(rep(c(0.1,0.5,1,2,3), each = 300),100),
         ncol = 300, byrow = TRUE),
  # early deactivation
  matrix(rep(exp(atan( ((300:1)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # late deactivation
  matrix(rep(exp(atan( ((300:1)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # early activation
  matrix(rep(exp(atan( ((1:300)-100)/50 )),50), ncol = 300, byrow = TRUE),
  # late activation
  matrix(rep(exp(atan( ((1:300)-200)/50 )),50), ncol = 300, byrow = TRUE),
  # transient
  matrix(rep(exp(atan( c((1:100)/33, rep(3,100), (100:1)/33) )),50), 
         ncol = 300, byrow = TRUE)
)
counts <- apply(means,2,function(cell_means){
  total <- rnbinom(1, mu = 7500, size = 4)
  rmultinom(1, total, cell_means)
})
rownames(counts) <- paste0('G',1:750)
colnames(counts) <- paste0('c',1:300)
sce <- SingleCellExperiment(assays = List(counts = counts))

data("slingshotExample")
rd <- slingshotExample$rd
cl <- slingshotExample$cl

dim(rd) # data representing cells in a reduced dimensional space
length(cl) # vector of cluster labels

# filter genes down to potential cell-type markers
# at least M (15) reads in at least N (15) cells
geneFilter <- apply(assays(sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
sce <- sce[geneFilter, ]

# quantile normalization
FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}
assays(sce)$norm <- FQnorm(assays(sce)$counts)

# reduce dimension by PCA
pca <- prcomp(t(log1p(assays(sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

# reduce dimension by UMAP
library(uwot)
rd2 <- uwot::umap(t(log1p(assays(sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

# add reduction results to SCE object
reducedDims(sce) <- SimpleList(PCA = rd1, UMAP = rd2)

# clustering by Gaussian mixture modeling
# BiocManager::install('mclust')
library(mclust)

cl1 <- Mclust(rd1)$classification
colData(sce)$GMM <- cl1

library(RColorBrewer)
plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)

# clustering by k-means
cl2 <- kmeans(rd1, centers = 4)$cluster # choose k=4 arbitarily
colData(sce)$kmeans <- cl2

plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

# run slingshot with the dimensionality reduction produced by PCA and cluster labels identified by Gaussian mixutre modeling
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')

# visualization: inferred lineage for the single-trajectory data with points colored by pseudotime
library(grDevices)
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

# see how the lineage structure was intially estimated by the cluster-based minimum spanning tree by using the type argument
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

# Identifying temporally dynamic genes
# BiocManager::install('tradeSeq')
library(tradeSeq)

# fit negative binomial GAM
sce <- fitGAM(sce)

# test for dynamic expression
ATres <- associationTest(sce)

# pick most dynamic 250 genes
topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])

# following analysis use example data
# identify global lineage structure

lin1 <- getLineages(rd, cl, start.clus = '1')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin1), lwd = 3, col = 'black')

# known endpoints can be specified
lin2 <- getLineages(rd, cl, start.clus= '1', end.clus = '3')

plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(lin2), lwd = 3, col = 'black', show.constraints = TRUE)

# constructing smooth curves
crv1 <- getCurves(lin1)
plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
lines(SlingshotDataSet(crv1), lwd = 3, col = 'black')

# For large datasets, we highgly recommend using the approx_points argument with slingshot (or getCurves)
# recommend a value of 100-200, hence the default value of 150
# restricting the number of unique points along the curve does not impose a similar limit on the number of unique pseudotime values, as demonstrated below
sce5 <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA',
                  approx_points = 5)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce5$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce5)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce5), lwd=2, col='black')

# construct multiple trajectory by adjusting omega