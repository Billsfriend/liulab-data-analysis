library(tidyverse)

read_lines('chr14.recal.vcf.gz', n_max = 40)[31:40]

korea1k <- 'ftp://biodisk.org/Release/Korea1K/SNV_Indel/Unrelated_BatchEffectFiltering/'

vcf_file <- read_delim(paste0(korea1k, 'chr14.recal.vcf.gz'),
                       comment = '#',
                       col_names = FALSE)

glimpse(vcf_file)

filter(vcf_file, X2 == 106204113) -> test

filter(vcf_file, X3 == 'rs1050501')

t <- filter(vcf_file, X3 == 'rs1050501')

t$X8

