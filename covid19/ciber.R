library(ggpubr)
library(tidyverse)

fraction <- read_tsv('../covid19/data/zhang-jiang2022inactivated/CIBERSORTx_Job2_Results.txt')

two_dose <- fraction %>%
  filter(str_detect(Mixture, '2nd')) %>%
  mutate(genotype = case_when(
    str_detect(Mixture, '2nd-2') ~ 'RR',
    TRUE ~ 'GG'
  )) %>%
  pivot_longer(cols = 2:7, names_to = 'cell_type', values_to = 'fraction_score')

two_dose %>%
  ggplot(aes(genotype, fraction_score)) +
  geom_jitter(width = 0.1) +
  stat_summary(geom = 'crossbar', fun = 'mean', width = 0.2, color = 'red') +
  facet_wrap(~ cell_type, scale = 'free') +
  font('x.text', size = 16)
