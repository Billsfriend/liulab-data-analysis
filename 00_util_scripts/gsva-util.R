library(GSVA)
library(GSVAdata)
library(Seurat)
library(pheatmap)
library(msigdbr)
library(limma)
library(tidyverse)
data(c2BroadSets)

# read in seurat dataset ---------
sobj <- read_rds('data/Wang2022/identImmune.rds') %>%
  SetIdent(value = "recessive")

cell_type_list <- SplitObject(sobj, split.by = 'cell_type')

# determine pathway database -----------
pathwaysDB <- msigdbr(species = "Homo sapiens",
               category = "C2",
               subcategory = 'CP:REACTOME') 
msigdbr_list <- split(x = pathwaysDB$gene_symbol,
                     f = pathwaysDB$gs_name)

# define util functions -----------
get_gsva_safely <- function(sub_sobj){
  try(get_gsva_result(sub_sobj))
}

get_gsva_result <- function(sub_sobj){
  name_of_element <- sub_sobj$SingleR.labels[1]
  sub_sobj %>%
    Idents() %>%
    factor(levels = c(FALSE, TRUE),
           ordered = F) ->
    grouping
  
  grouping_mtx <- model.matrix(~grouping)
  
  colnames(grouping_mtx) <- levels(grouping)
  
  sub_sobj %>%
    GetAssayData() %>%
    gsva(msigdbr_list,
         parallel.sz = 10,
         min.sz = 10,
         max.sz = 100) %>%
    lmFit(grouping_mtx) %>%
    eBayes() %>%
    topTable(adjust = 'fdr',
             coef = 2,
             number = Inf,
             resort.by = 't',
             p.value = 0.05)
}

# test for element-1 in list --------
gsva_result_1 <- get_gsva_result(cell_type_list[[1]])

gsva_result_list <- lapply(cell_type_list,
                           get_gsva_safely)

# write result to disk
write.csv(gsva_result_1, 'results/GSVA/wang2022/epi_gsva.csv')

for (i in 1:6) {
  write.csv(gsva_result_list[[i]],
            paste0('results/GSVA/wang2022/immune/',
                   names(gsva_result_list)[i],
                   '.csv'))
}

# plot digested pathways for GSVA ------
gsva_result_1 %>%
  rownames_to_column('pathway') %>%
  mutate(pathway = str_to_sentence(pathway)) %>%
  top_n(7, t) %>%
  ggplot(aes(x = pathway, y = t)) +
  geom_col(aes(fill = adj.P.Val)) +
  coord_flip()


# draw barplot of GSVA t value
ggplot(richpath, aes(pathway, t))+
  geom_col(aes(fill=factor)) +
  theme_classic()+
  theme(
    text = element_text(size = 16), 
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank()
  )+
  geom_text(aes(label=pathway, hjust=0.5+0.5*t/abs(t)), size=5,nudge_y = -1*(richpath$t)-0.1*richpath$t/abs(richpath$t))+
  coord_flip() +
  geom_hline(yintercept=0)+
  ylab('t value of GSVA analysis')+
  ggtitle(richpath$subset[1])

# proportion of T cell subset in all T cells
natperc <- read_delim('results/GSVA/digest/naturepercent.txt')

natperc$subset <- reorder(natperc$subset, natperc$mean, mean)

nperc <- ggplot(natperc,aes(subset,mean))
nperc+
  geom_col(aes(fill=factor))+
  coord_flip()

# generate heatmap by geom_tile
ggplot(natperc, aes(x=factor,y=subset,fill=mean))+
  geom_tile()+
  geom_raster()+
  scale_fill_viridis_c()+
  labs(fill='Proportion in all T cells')+
  theme(text=element_text(size=16))
