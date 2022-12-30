# 2022.3.23
# re-analyze smart-seq TPM data with Seurat
# from 2020 Zhang Cell paper 

# edit in 2022.9.28

library(data.table)
library(Seurat)
library(tidyverse)

# load smart-seq TPM data ---------
data <- fread('data/cell/CRC.Leukocyte.Smart-seq2.TPM.txt.gz')
data[1:5,1:5]

genes <- data$V1

genotyped_smart <- c('P0104','P0411','P1212','P0305','P0720','P0825')
II_smart <- c('P0104','P0411','P1212')
IT_smart <- c('P0305','P0720','P0825')

select(data, contains(genotyped_smart)) -> data

sobj <- CreateSeuratObject(data,
                           row.names = genes,
                           names.field = 3)

sobj@meta.data -> seurat_meta

seurat_meta$CellName <- rownames(seurat_meta)

# load metadata
smrt_meta <- fread("data/cell/CRC.Leukocyte.Smart-seq2.Metadata.txt")

# assign I232T genotype in smrt_meta
smrt_meta %>%
  mutate(genotype = case_when(
    Sample %in% II_smart ~ 'II',
    Sample %in% IT_smart ~ 'IT',
    TRUE ~ 'NA'
  )) %>%
  filter(genotype != 'NA') %>%
  mutate(rowname = CellName) %>%
  as.data.frame() %>%
  right_join(seurat_meta) %>%
  column_to_rownames('rowname') -> smrt_meta

sobj@meta.data <- smrt_meta

# save genotyped data
write_rds(sobj, 'data/cell/tpmSmart.rds')

# Visualize the number of cell counts per sample
sobj@meta.data %>% 
  ggplot(aes(x=genotype, fill=genotype)) + 
  geom_bar() +
  theme_classic() +
  labs(title = 'Number of cells')

# Visualize the number UMIs/transcripts per cell
sobj@meta.data %>% 
  ggplot(aes(color=genotype, x=filter.nUMI, fill=genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# load 10x TPM data ------------
data <- fread("CRC-I/data/Zhang-Yu-2020/CRC.Leukocyte.10x.TPM.txt.gz")

data[1:5, 1:5]

genes <- data$V1

II_10x <- c('P0410','P0323','P0408','P1026','P0104')
IT_10x <- c('P0123','P0202','P0613','P1025','P0305')

# the cell naming: P_N_P0410_01391, first letter is always 'P', second letter is N/P/T representing tissue source.
# here we only need cells in tumor
select(data, contains(c(II_10x, IT_10x))) -> data
select(data, contains('_T_')) -> data

sobj <- CreateSeuratObject(data,
                           row.names = genes,
                           names.field = 3,
                           min.cells = 1,
                           min.features = 1)

sobj@meta.data -> seurat_meta

seurat_meta$CellName <- rownames(seurat_meta)

# load metadata
tenx_meta <- read_delim("CRC-I/data/Zhang-Yu-2020/CRC.Leukocyte.10x.Metadata.txt")

# assign I232T genotype in tenx_meta
tenx_meta %>%
  mutate(genotype = case_when(
    Sample %in% II_10x ~ 'II',
    Sample %in% IT_10x ~ 'IT',
    TRUE ~ 'NA'
  )) %>%
  filter(genotype != 'NA') %>%
  right_join(seurat_meta) %>%
  column_to_rownames('CellName') -> tenx_meta

sobj@meta.data <- tenx_meta
 
# save 10x seurat file
write_rds(sobj, 'CRC-I/data/Zhang-Yu-2020/zy2020_tumor10x.rds')

# Visualize the number of cell counts per sample
sobj@meta.data %>% 
  ggplot(aes(x=genotype, fill=genotype)) + 
  geom_bar() +
  theme_classic() +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("Number of cells")

# Visualize the number UMIs/transcripts per cell
sobj@meta.data %>% 
  ggplot(aes(color=genotype, x=filter.nUMI, fill=genotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)

# Visualize the distribution of genes detected per cell via histogram
sobj@meta.data %>% 
  ggplot(aes(color=genotype, x=filter.nGene, fill= genotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
sobj@meta.data %>% 
  ggplot(aes(x=genotype, y=log10(filter.nGene), fill=genotype)) + 
  geom_violin() + 
  theme_classic() +
  ggtitle("Number of genes per cell")+
  BoldTitle()

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
sobj@meta.data %>% 
  ggplot(aes(x=filter.nUMI, y=filter.nGene)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250)

# Visualize the overall ribo.per
sobj@meta.data %>%
  ggplot(aes(x=ribo.per, color = genotype, fill=genotype)) +
  geom_density(alpha = 0.2) +
  theme_classic()
