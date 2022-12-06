# topGO vignette ver 2021 Oct

# BiocManager::install('topGO')
# BiocManager::install('ALL')
# BiocManager::install('hgu95av2.db')
library(topGO)
library(ALL) # leukemia microarray data in Chiaretti 2004
data(ALL) # expressionSet object
data(geneList) # p values of DEGs

affyLib <- paste(annotation(ALL), "db", sep = ".") # name of microarray used
library(package = affyLib, character.only = TRUE)

sum(topDiffGenes(geneList)) # topDiffGenes() only return p value which < 0.01

sampleGOdata <- new('topGOdata', description = 'Simple session', ontology = 'BP', allGenes = geneList, geneSel = topDiffGenes, nodeSize = 10, annot = annFUN.db, affyLib = affyLib)
# nodeSize = 10 means GO terms with less than 10 annotated genes will be pruned out

sampleGOdata # see summary of data

resultFisher <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'fisher')
# classic GO algorithm each GO category is tested independently
# fisher exact test is based on gene counts

resultFisher # see summary of Fisher test

resultKS <- runTest(sampleGOdata, algorithm = 'classic', statistic = 'KS')
# Kolmogorov-Smirnov test compute score for each genes

resultKS.elim <- runTest(sampleGOdata, algorithm = 'elim', statistic = 'KS')
# elim algorithm is used here
# topGO use weight01 as default algorithm

allRes <- GenTable(sampleGOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, orderBy = 'elimKS', ranksOf = 'classicFisher', topNodes = 10)
# generate a result table, colnames and coldata need to be specified manually




