# https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

# BiocManager::install('fgsea')

library(fgsea)
library(data.table)
library(ggplot2)

data("examplePathways") # mouse Reactome
data("exampleRanks") # mouse Th1 polarization
set.seed(42)

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500)

head(fgseaRes[order(pval), ])

fgseaRes <- fgsea(pathways = examplePathways, 
                  stats    = exampleRanks,
                  minSize  = 15,
                  maxSize  = 500,
                  eps = 0.0) # set eps = 0.0 so that p value is more accurate

# make an enrichment plot of a pathway
plotEnrichment(examplePathways[["5991130_Programmed_Cell_Death"]],
               exampleRanks) + labs(title="Programmed Cell Death")

# or a table plot for a selected bunch of pathways
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
plotGseaTable(examplePathways[topPathways], exampleRanks, fgseaRes, 
              gseaParam=0.5)

# to collapse very similar pathways into one entry
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                      examplePathways, exampleRanks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
  order(-NES), pathway]
plotGseaTable(examplePathways[mainPathways], exampleRanks, fgseaRes, 
              gseaParam = 0.5)

# save the result in text
fwrite(fgseaRes, file="results/fgseaRes.txt", sep="\t", sep2=c("", " ", ""))

# make results more readable by annotate gene symbols
library(org.Mm.eg.db)
fgseaResMain <- fgseaRes[match(mainPathways, pathway)]
fgseaResMain[, leadingEdge := mapIdsList(
  x=org.Mm.eg.db, 
  keys=leadingEdge,
  keytype="ENTREZID", 
  column="SYMBOL")]
fwrite(fgseaResMain, file="fgseaResMain.txt", sep="\t", sep2=c("", " ", ""))

# choose reactome pathways from reactome.db package
# BiocManager::install('reactome.db')
pathways <- reactomePathways(names(exampleRanks))
fgseaRes <- fgsea(pathways, exampleRanks, maxSize=500)
head(fgseaRes[order(pval), ])

# or, start from files
rnk.file <- system.file("extdata", "naive.vs.th1.rnk", package="fgsea")
gmt.file <- system.file("extdata", "mouse.reactome.gmt", package="fgsea")

# loading ranks
ranks <- read.table(rnk.file,
                    header=TRUE, colClasses = c("character", "numeric"))
ranks <- setNames(ranks$t, ranks$ID)
str(ranks)

# loading pathways
pathways <- gmtPathways(gmt.file)
str(head(pathways))
