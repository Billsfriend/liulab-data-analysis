# haven: import data from different data base file type
library(haven)
library(data.table)
crc <- read_sav('data/1006crc.sav')
fread('data/1006CRC-SNP.sav')
