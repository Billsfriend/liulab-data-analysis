# GSVA package vignette
# https://www.bioconductor.org/packages/release/bioc/vignettes/GSVA/inst/doc/GSVA.html

BiocManager::install('GSVA')
BiocManager::install('GSVAdata')
library(GSVA)

p <- 10000 ## number of genes
n <- 30    ## number of samples
## simulate expression values from a standard Gaussian distribution
X <- matrix(rnorm(p*n), nrow=p,
            dimnames=list(paste0("g", 1:p), paste0("s", 1:n)))
X[1:5, 1:5]

## sample gene set sizes
gs <- as.list(sample(10:100, size=100, replace=TRUE))
## generate sample gene sets
gs <- lapply(gs, function(n, p)
  paste0("g", sample(1:p, size=n, replace=FALSE)), p)
names(gs) <- paste0("gs", 1:length(gs))

# example of how to perform GSVA on expression matrix and gene sets
gsva.es <- gsva(X, gs, verbose=FALSE)
dim(gsva.es)

gsva.es[1:5, 1:5]

# load library for standardized gene set object class
library(GSEABase)
# library for human gene set database
library(GSVAdata)

data(c2BroadSets)
class(c2BroadSets)

# load example published data to analyze
library(Biobase)

data(commonPickrellHuang)

stopifnot(identical(sampleNames(huangArrayRMAnoBatchCommon_eset),
                    sampleNames(pickrellCountsArgonneCQNcommon_eset)))

# pick gene sets we want
canonicalC2BroadSets <- c2BroadSets[c(grep("^KEGG", names(c2BroadSets)),
                                      grep("^REACTOME", names(c2BroadSets)),
                                      grep("^BIOCARTA", names(c2BroadSets)))]
data(genderGenesEntrez)

MSY <- GeneSet(msYgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="MSY")
XiE <- GeneSet(XiEgenesEntrez, geneIdType=EntrezIdentifier(),
               collectionType=BroadCollection(category="c2"), setName="XiE")
canonicalC2BroadSets <- GeneSetCollection(c(canonicalC2BroadSets, MSY, XiE))

# perform GSVA for published data
# for integer data kcdf is 'Gaussian' for default
esmicro <- gsva(huangArrayRMAnoBatchCommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500)
# for TPM-like data (not log-TPM), kcdf = 'Poisson'
esrnaseq <- gsva(pickrellCountsArgonneCQNcommon_eset, canonicalC2BroadSets, min.sz=5, max.sz=500,
                 kcdf="Poisson")


library(edgeR)

# perform log transformation on Counts per million data
lcpms <- cpm(exprs(pickrellCountsArgonneCQNcommon_eset), log=TRUE)

# calculate Spearman correlations between gene expression profiles
# of the previous log-CPM values and the microarray RMA values
genecorrs <- sapply(1:nrow(lcpms),
                    function(i, expmicro, exprnaseq) cor(expmicro[i, ], exprnaseq[i, ],
                                                         method="spearman"),
                    exprs(huangArrayRMAnoBatchCommon_eset), lcpms)
names(genecorrs) <- rownames(lcpms)

# calculate Spearman correlations between GSVA enrichment scores 
# derived from the microarray and the RNA-seq data
pwycorrs <- sapply(1:nrow(esmicro),
                   function(i, esmicro, esrnaseq) cor(esmicro[i, ], esrnaseq[i, ],
                                                      method="spearman"),
                   exprs(esmicro), exprs(esrnaseq))
names(pwycorrs) <- rownames(esmicro)

# show histograms of these correlations
par(mfrow=c(1, 2), mar=c(4, 5, 3, 2))
hist(genecorrs, xlab="Spearman correlation", main="Gene level\n(RNA-seq log-CPMs vs microarray RMA)",
     xlim=c(-1, 1), col="grey", las=1)
hist(pwycorrs, xlab="Spearman correlation", main="Pathway level\n(GSVA enrichment scores)",
     xlim=c(-1, 1), col="grey", las=1)

# differential expression analysis of a published leukemia data
data(leukemia)

# calculate GSVA enrichment scores setting the minimum and maximum 
# gene set size to 10 and 500 genes, respectively
leukemia_es <- gsva(leukemia_eset, c2BroadSets, min.sz=10, max.sz=500)

# use commands from limma package to analyze between leukemia subtypes
mod <- model.matrix(~ factor(leukemia_es$subtype))
colnames(mod) <- c("ALL", "MLLvsALL")
fit <- lmFit(leukemia_es, mod)
fit <- eBayes(fit)
res <- decideTests(fit, p.value=0.01)
summary(res)

# draw volcano plot of 96 DEG pathways
tt <- topTable(fit, coef=2, n=Inf)
DEpwys <- rownames(tt)[tt$adj.P.Val <= 0.01]
plot(tt$logFC, -log10(tt$P.Value), pch=".", cex=4, col=grey(0.75),
     main="", xlab="GSVA enrichment score difference", ylab=expression(-log[10]~~Raw~P-value))
abline(h=-log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), col=grey(0.5), lwd=1, lty=2)
points(tt$logFC[match(DEpwys, rownames(tt))],
       -log10(tt$P.Value[match(DEpwys, rownames(tt))]), pch=".", cex=5, col="darkred")
text(max(tt$logFC)*0.85, -log10(max(tt$P.Value[tt$adj.P.Val <= 0.01])), "1% FDR", pos=3)

# heatmap of 96 DEG pathways
DEpwys_es <- exprs(leukemia_es[DEpwys, ])
colorLegend <- c("darkred", "darkblue")
names(colorLegend) <- c("ALL", "MLL")
sample.color.map <- colorLegend[pData(leukemia_es)[, "subtype"]]
names(sample.color.map) <- colnames(DEpwys_es)
sampleClustering <- hclust(as.dist(1-cor(DEpwys_es, method="spearman")),
                           method="complete")
geneSetClustering <- hclust(as.dist(1-cor(t(DEpwys_es), method="pearson")),
                            method="complete")
heatmap(DEpwys_es, ColSideColors=sample.color.map, xlab="samples",
        ylab="Pathways", margins=c(2, 20),
        labRow=substr(gsub("_", " ", gsub("^KEGG_", "",
                                          rownames(DEpwys_es))), 1, 35),
        labCol="", scale="row", Colv=as.dendrogram(sampleClustering),
        Rowv=as.dendrogram(geneSetClustering))
legend("topleft", names(colorLegend), fill=colorLegend, inset=0.01, bg="white")

m_df = msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") #选取物种人类
