library(Seurat)
library(Platypus)
library(pheatmap)
library(GSVA)
library(GSVAdata)
library(msigdbr)
library(limma)
library(tidyverse)

# read in seurat data ----------
sobj <- read_rds('data/zhang_T_sample.rds')

sobj$recessive <- FALSE

sobj$recessive[str_detect(sobj$genotype, 'TT')] <- TRUE

sobj <- SetIdent(sobj, value = 'recessive')

write_rds(sobj, 'data/zhang_B_cell.rds')

SplitObject(sobj, split.by = 'celltype') -> celltype_list

celltype_list$Plasma %>%
  DotPlot(features = c('IGHG1',
                       'IGHG2',
                       'IGHG3',
                       'IGHG4',
                       'IGHA1',
                       'IGHA2',
                       'IGHE',
                       'IGHM'),
          group.by = 'genotype')+
  labs(x = 'genes', y = 'genotype')

findmarkInTT <- function(srt){try(FindMarkers(
  srt,
  group.by = 'genotype',
  ident.1 = 'TT',
  fc.name = "avg_logFC"))
}

degList <- map(celltype_list, findmarkInTT)

degList$TEFF %>%
  GEX_volcano(input.type = 'findmarkers',
              condition.1 = 'TT',
              condition.2 = 'II+IT',
              n.label.up = 6)+
  labs(subtitle = 'effector T cell')

write_rds(degList, 'results/deg_zhang_myeloid.rds')

# gsva -------

pathwaysDB <- msigdbr(species = "Homo sapiens",
                      category = "C2",
                      subcategory = 'CP:REACTOME')

# test for element-1 in list --------
gsva_result_1 <- get_gsva_result(celltype_list[[1]])

gsva_result_list <- lapply(celltype_list,
                           try(get_gsva_result))

# visualization of digest GSVA --------
dgst <- read_tsv('results/zhang_gsva_digest.txt')

dgst %>%
  filter(cell == 'tDC' & t > 0) %>%
  ggplot()+
  geom_col(aes(x = pathway, y = t, fill = adj.P.Val))+
  coord_flip()+
  scale_fill_distiller(palette = 'Oranges')+
  theme(text = element_text(size = 18))+
  labs(y = 't value of GSVA',
       title = 'TT tolerogenic DC cells')

# downregulated pathways
dgst %>%
  filter(cell == 'tDC' & t < 0) %>%
  ggplot()+
  geom_col(aes(x = pathway, y = t, fill = adj.P.Val))+
  coord_flip()+
  scale_fill_distiller(palette = 'Blues')+
  theme(text = element_text(size = 14))+
  labs(y = 't value of GSVA',
       title = 'TT tolerogenic DCs')
