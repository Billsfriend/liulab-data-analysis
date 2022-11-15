library(biodbUniprot)
library(biodb)
library(tidyverse)

mybiodb <- newInst()

conn <- mybiodb$getFactory()$createConn('uniprot')

entries <- conn$getEntry(c('P01011', 'P09237'))

mybiodb$entriesToDataframe(entries)

# BiocManager::install('UniProt.ws')
library(UniProt.ws)

up <- UniProt.ws()
