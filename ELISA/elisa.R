library(tidyverse)
library(ggpubr)

# data 221105 --------------
read_delim('ELISA/data/elisa-fcgr2b-ova-z93-d14.txt') -> data

data %>%
  filter(subtype %in% c('IgG1','IgG2c')) %>%
  ggplot(aes(group, y = titer, color = group))+
  stat_summary(geom = 'col', fun = 'mean', fill = 'white', linewidth = 1)+
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.3)+
  geom_jitter(width = 0.1, height = 0, size = 2)+
  expand_limits(y = 0)+
  ggprism::theme_prism()+
  stat_compare_means(comparisons = list(c('II','IT')), method = 't.test')+
  scale_color_discrete(label = c('II (n=6)','IT (n=5)'))+
  facet_wrap(~subtype)+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))
  
# data 221226 -----------
read_delim('ELISA/data/elisa-fcgr2b-hel-aln-d7-17.txt') -> data

data %>%
  pivot_longer(3:4, names_to = 'group', values_to = 'titer') %>%
  filter(!is.na(titer)) %>%
  mutate(day = str_glue('d{day}')) %>%
  filter(day == 'd17') %>%
  ggplot(aes(group, titer, color = group)) +
  stat_summary(geom = 'col', fun = 'mean', fill = 'white', linewidth = 1)+
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 0.3)+
  geom_jitter(width = 0.1, height = 0, size = 2)+
  expand_limits(y = 0)+
  ggprism::theme_prism()+
  stat_compare_means(comparisons = list(c('II','IT')), method = 't.test')+
  scale_color_discrete(label = c('II (n=6)','IT (n=5)'))+
  facet_wrap(~type+day, scales = 'free')+
  scale_y_continuous(expand = expansion(mult = c(0,0.1)))

