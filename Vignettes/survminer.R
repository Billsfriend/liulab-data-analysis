library(survminer)
library(survival)
library(tidyverse)

attach(lung) # attach example data from 'survival' package
Surv(time, status) # create Surv object

fit <- survfit(Surv(time, status) ~ sex, # surv time and status, and grouping factor
  data = lung
)

ggsurvplot(fit,
  data = lung,
  pval = TRUE, # add p value
  risk.table = TRUE, # add risk table
  surv.median.line = "hv", # add median line
  conf.int = TRUE
) # add confidence interval

ggsurvplot(fit,
  data = lung,
  conf.int = TRUE,
  fun = "cumhaz"
) # plot cumulative hazard curve

# this function will add curve of group 'all samples'
ggsurvplot_add_all(fit,
  data = lung,
  palette = "rickandmorty"
) # palette can be changed

ggsurvplot(fit,
  data = lung,
  xlab = "Follow up time (d)", # x axis label
  legend.title = "", # remove legend title
  legend.labs = c("Male", "Female"), # set legend text
  break.x.by = 100, # set tick break on x axis
  pval = TRUE,
  risk.table = TRUE,
  surv.median.line = "hv",
  conf.int = TRUE
)


# with my data of CRC -----------
crc <- read_delim("data/other/crcSurvival-cross.txt")

crc %>%
  mutate(group = case_when(
    I232T != 3 & G396R == 3 ~ 'survival++',
    I232T == 3 & G396R != 3 ~ 'survival--',
    TRUE ~ 'survival+-'
  )) -> crc

clinVar <- read_csv('data/other/1006CRC-20201120-AllVar.csv')

clinVar %>%
  mutate(patient = str_remove(patient ,'\\.')) %>%
  filter(survival_month != '#NULL!') %>%
  left_join(crc, by = c('patient' = 'id')) ->
  clinMerge

survfit(Surv(crc$month, crc$OS) ~ crc$group) %>%
  ggsurvplot(data = crc,
             xlab = "Follow up time (month)",
             legend.labs = c("bad + bad",
                             "good + bad",
                             "good + good"),
             pval = TRUE,
             pval.size = 12,
             break.x.by = 20,
             risk.table = TRUE,
             surv.median.line = "hv",
             title = "CRC patients OS with FCGR2B-IT and IGHG1-GR combination",
             font.legend = 14)

crcos <- survfit(Surv(clinVar$month, clinVar$OS) ~ clinVar$G396R)

ggsurvplot(crcos,
  data = clinVar,
  xlab = "Follow up time (month)",
  #legend.labs = c("II", "IT", "TT"),
  pval = TRUE,
  pval.size = 12,
  break.x.by = 20,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients OS with FCGR2b-I232T"
)

ggsurvplot(crcos,
  data = crc,
  fun = "cumhaz",
  xlab = "Follow up time (month)",
  legend.labs = c("II", "IT", "TT"),
  pval = TRUE,
  conf.int = TRUE
)

crcpfs <- survfit(Surv(crc$month, crc$PFS) ~ crc$G396R)

ggsurvplot(crcpfs,
  data = crc,
  xlab = "Follow up time (month)",
  legend.labs = c("II", "IT", "TT"),
  pval = TRUE,
  pval.size = 12,
  break.x.by = 20,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients PFS with FCGR2b-I232T"
)

# evaluate combination with IGHG1-G396R
library(data.table)
setDT(crc)

crc$goodgene <- 2
crc[I232T < 3 & G396R == 3, goodgene := 1] # create a new column based on existed values
crc[I232T == 3 & G396R < 3, goodgene := 3]


crcpfsCross <- survfit(Surv(crc$month, crc$PFS) ~ crc$goodgene)

ggsurvplot(crcpfsCross,
  data = crc,
  xlab = "Follow up time (month)",
  legend.labs = c("Survival+/+", "Survival+/-", "Survival-/-"),
  pval = TRUE,
  pval.size = 12,
  break.x.by = 20,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients PFS with FCGR2B-I232T"
)

crcosCross <- survfit(Surv(crc$month, crc$OS) ~ crc$goodgene)

ggsurvplot(crcosCross,
  data = crc,
  xlab = "Follow up time (month)",
  legend.labs = c("Survival++", "Survival+", "Survival-"),
  pval = TRUE,
  pval.size = 12,
  break.x.by = 20,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients OS with FCGR2B-I232T"
)

crc[I232T < 3 & G396R == 3, good2gene := 1]
crc[I232T < 3 & G396R < 3, good2gene := 2] # create a new column based on existed values
crc[I232T == 3 & G396R == 3, good2gene := 3]
crc[I232T == 3 & G396R < 3, good2gene := 4]

crcosCross2 <- survfit(Surv(crc$month, crc$OS) ~ crc$good2gene)

ggsurvplot(crcosCross2,
  data = crc,
  xlab = "Follow up time (month)",
  legend.labs = c("FCGR2B+ IGHG1+",
                  "FCGR2B+ IGHG1-",
                  "FCGR2B- IGHG1+",
                  "FCGR2B- IGHG1-"),
  palette = c('blue', 'green', 'orange', 'red'),
  pval = TRUE,
  pval.size = 12,
  break.x.by = 20,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients OS with FCGR2B-I232T"
)

# combine with MSI variable ----------
clinMerge %>%
  filter((MSI == '0' & SNP1 == 3)|(MSI == '1' & SNP1 == 1) ) ->
  crc

crc_OS_MSI <- survfit(Surv(crc$survival_month, crc$OS.x) ~ crc$MSI + crc$SNP1)

ggsurvplot(crc_OS_MSI,
           data = crc,
           xlab = "Follow up time (month)",
             legend.labs = c("RR + MSI",
                             "GG + MSS"),
           # #palette = c('blue', 'green', 'orange', 'red'),
           pval = TRUE,
           pval.size = 12,
           break.x.by = 20,
           risk.table = 'nrisk_cumcensor',
           title = "CRC patients OS",
           tables.height = 0.35
)

# with TCGA dataset -----------
coadSurvival <- read_csv("results/coad+read_survival.csv")

coadOS <- survfit(Surv(crcSurvCluster$days_to_last_follow_up, crcSurvCluster$dead) ~ crcSurvCluster$CD4_CTLA4)

ggsurvplot(coadOS,
  data = crcSurvCluster,
  xlab = "Follow up time (days)",
  # legend.labs = c('II','IT','TT'),
  pval = TRUE,
  pval.size = 12,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "CRC patients OS with FCGR2b-I232T"
)

ggsurvplot(coadOS,
           data = crcSurvival)
