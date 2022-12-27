# BiocManager::install('ChIPseeker')
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(clusterProfiler)
library(tidyverse)

files <- getSampleFiles()
peak <- readPeakFile(files[[4]])
peak

# coverage plot
covplot(peak, weightCol="V5")

covplot(peak,
        weightCol="V5",
        chrs = c('chr1'),
        xlim = c(634024, 634025))

capSummit <- readPeakFile('../data/capsid_summits.bed')
inpSummit <- readPeakFile('../data/input_summits.bed')

covplot(capSummit,
        weightCol = 'V5',
        chrs = c('chr1'))

covplot(inpSummit,
        weightCol = 'V5',
        chrs = c('chr1'))


# get a matrix of promoter regions
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb = txdb,
                         upstream = 3000,
                         downstream = 3000)

tagMatrix <- getTagMatrix(peak,
                          windows = promoter)

capsidPromoter <- getTagMatrix(capSummit,
                               windows = promoter)

inputPromoter <- getTagMatrix(inpSummit,
                               windows = promoter)

# heatmap of binding to promoter
tagHeatmap(tagMatrix,
           xlim=c(-3000, 3000))

tagHeatmap(capsidPromoter,
           xlim=c(-3000, 3000))

# plot heatmap directly from bed file
peakHeatmap(files[[4]],
            TxDb=txdb,
            upstream=3000,
            downstream=3000, color="red")

# Average profile of peaks binding to promoter
plotAvgProf(tagMatrix,
            xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency")

plotAvgProf(capsidPromoter,
            xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency")

plotAvgProf(inputPromoter,
            xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')",
            ylab = "Read Count Frequency")

# add confidence interval
plotAvgProf(tagMatrix,
            xlim=c(-3000, 3000),
            conf = 0.95,
            resample = 1000)

## The results of binning method and normal method are nearly the same.
## get TSS by gene
tagMatrix_binning <- getTagMatrix(
  peak = peak,
  TxDb = txdb, 
  upstream = 3000,
  downstream = 3000, 
  type = "start_site",
  by = "gene", 
  weightCol = "V5",
  nbin = 800)

## Here uses `plotPeakProf2` to do all things in one step.
## Gene body regions having lengths smaller than nbin will be filtered
## A message will be given to warning users about that.
## >> 9 peaks(0.872093%), having lengths smaller than 800bp, are filtered...

## the ignore_strand is FALSE in default. We put here to emphasize that.
## We will not show it again in the below example
plotPeakProf2(peak = peak,
              upstream = rel(0.2),
              downstream = rel(0.2),
              conf = 0.95,
              by = "gene",
              type = "body",
              nbin = 800,
              TxDb = txdb,
              weightCol = "V5",
              ignore_strand = F)

# get the profile ChIP peaks binding to gene body regions with no flank extension or flank extension decided by actual length.
## The first method using getBioRegion(), getTagMatrix() and plotPeakProf() to plot in three steps.
genebody <- getBioRegion(TxDb = txdb,
                         by = "gene",
                         type = "body")

matrix_no_flankextension <- getTagMatrix(
  peak,
  windows = genebody,
  nbin = 800)

plotPeakProf(matrix_no_flankextension,
             conf = 0.95)

## The second method of using getTagMatrix() and plotPeakProf() to plot in two steps
matrix_actual_extension <- getTagMatrix(
  peak,
  windows = genebody,
  nbin = 800,
  upstream = 1000,
  downstream = 1000)

plotPeakProf(matrix_actual_extension,
             conf = 0.95)

# peak annotation ----------
peakAnno <- annotatePeak(capSummit,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")

inpAnno <- annotatePeak(inpSummit,
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb,
                         annoDb = "org.Hs.eg.db")

# visualize peak annotation
plotAnnoPie(peakAnno)
plotAnnoPie(inpAnno)

plotAnnoBar(peakAnno)
plotAnnoBar(inpAnno)

vennpie(peakAnno)

upsetplot(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)

# binding sites distances to TSS
plotDistToTSS(peakAnno,
              title = "Distribution of transcription factor-binding loci\nrelative to TSS")

# seq2gene to find nearest genes from range data
gene <- seq2gene(peak,
                 tssRegion = c(-1000, 1000),
                 flankDistance = 3000,
                 TxDb = txdb)

# peak data comparison -----------
tagMatrixList <- lapply(files,
                        getTagMatrix,
                        windows = promoter)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
plotAvgProf(tagMatrixList,
            xlim=c(-3000, 3000),
            conf=0.95,
            resample=500,
            facet="row")

peakAnnoList <- lapply(files,
                       annotatePeak,
                       TxDb=txdb,
                       tssRegion=c(-3000, 3000),
                       verbose=FALSE)

# venn plot of binding genes
genes = lapply(c(peakAnno, inpAnno),
              function(i) as.data.frame(i)$geneId)
vennplot(genes)
