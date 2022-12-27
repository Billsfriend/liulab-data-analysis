library(tidyverse)
library(ggpubr)
library(survminer)
library(survival)

read_tsv('CRC-I/data/other/salmonella-mice-220926.txt') %>%
  select(-id)-> data

read_tsv('CRC-I/data/other/jax-weight.txt') %>%
  add_column(group = 'jax')-> jax

# calculate proportion change --------
# absolute weight at beginning
data %>%
  ggplot(aes(x = group, y = `0`)) +
  geom_jitter(width = 0.01, height = 0) +
  stat_summary(geom = 'crossbar', fun = 'mean', color = 'red', width = 0.2) +
  facet_wrap(~sex) +
  labs(title = 'mice weight at day 0', y = 'weight (g)') +
  theme_pubr() +
  stat_compare_means()

data %>%
  mutate(across(where(is.numeric), ~ . / `0`)) %>%
  pivot_longer(cols = where(is.numeric), names_to = 'time', values_to = 'proportion') %>%
  filter(!is.na(proportion)) %>%
  mutate(time = as.numeric(time)) -> data_prop

data %>%
  pivot_longer(cols = where(is.numeric), names_to = 'time', values_to = 'weight') %>%
  filter(!is.na(weight)) %>%
  mutate(time = as.numeric(time))-> data_abs

data_abs %>%
  group_by(group, sex, time) %>%
  summarise(mean = mean(weight), sd = sd(weight)) %>%
  bind_rows(jax)->
  data_mean_sd

# prism-like time-data graph
data_prop %>%
  filter(sex == 'Female') %>%
  ggplot(aes(x = time, y = proportion, group = group, color = group))+
  stat_summary(geom = 'point', fun = 'mean', size = 2) +
  stat_summary(geom = 'errorbar', fun.data  = 'mean_sd', width = 1, linewidth = 1) +
  stat_summary(geom = 'path', fun = 'mean', linewidth = 1) +
  geom_hline(yintercept = 1, linetype = 'dashed')+
  geom_vline(xintercept = 41, linetype = 'dotted')+
  labs(title = 'Mice weight change (female)',
       y = 'proportion to original weight',
       x = 'time (days)') +
  ggprism::theme_prism() +
  theme(legend.position = 'top')+
  stat_compare_means(label = 'p.signif')+
  scale_color_discrete(label = c('IgG1 GR (n=4)','WT (n=5)'))

# plot absolute weight
data_abs %>%
  filter(sex == 'Female') %>%
  ggplot(aes(x = time, y = weight, color = group))+
  stat_summary(geom = 'point', fun = 'mean', size = 2) +
  stat_summary(geom = 'errorbar', fun.data  = 'mean_se', width = 1, size = 1) +
  stat_summary(geom = 'path', fun = 'mean', size = 1) +
  stat_compare_means(label = 'p.signif')+
  labs(title = 'Mice weight change (male)',
       x = 'time (days)') +
  ggprism::theme_prism() +
  theme(legend.position = 'top')+
  scale_color_discrete(label = c('IgG1 GR (n=5)','WT (n=6)'))

compare_means(weight ~ group ,data_abs, group.by = c('sex', 'time')) %>%
  filter(p < 0.05)

data_mean_sd %>%
  filter(sex == 'Female' & time < 50) %>%
  ggplot(aes(time, mean, color = group))+
  geom_path(linewidth = 1)+
  geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), linewidth = 1)+
  #facet_wrap()+
  ggprism::theme_prism()+
  labs(x = 'days post infection',
       y = 'weight (g)',
       title = 'Female mice after 1e5 Salmonella infection')

# survival anaylsis --------
suv <- read_delim('data/other/salmonella-survival.txt')

# use surv_summary + ggsurvplot_df to more flexibly plot survival curve
suv %>%
  filter(!str_detect(group, 'IgG2c')) ->
  suv_set

survfit(Surv(suv_set$days, suv_set$dead) ~ suv_set$group) %>%
  ggsurvplot(
    data = suv_set,
    pval = TRUE,
    xlab = "Follow up time (days)",
    legend.labs = c("IgG1-GR (n=13)",
                    "IgG1-tailless (n=8)",
                    "IgG1-TMKO (n=14)",
                    'WT (n=12)'),
    risk.table = TRUE)

suv %>%
  filter(!str_detect(group, 'IgG1')) ->
  suv_set2

survfit(Surv(suv_set2$days, suv_set2$dead) ~ suv_set2$group) %>%
  ggsurvplot(
    data = suv_set2,
    pval = TRUE,
    xlab = "Follow up time (days)",
    legend.labs = c("IgG2c-GR (n=4)",
                    "IgG2c-tailless (n=13)",
                    'WT (n=12)'))



# for ELISA results ------------
filter(total_length, group != 'PC') -> total_length

ggplot(total_length, aes(group, IGM, color = group)) +
  geom_jitter(width = 0.05, height = 0, size = 3) +
  stat_summary(geom = 'crossbar', width = 0.3, fun = 'mean', color = 'red') +
  theme_pubr(legend = 'none') +
  labs(title = 'OVA-specific IgM level', y = 'titer') +
  ggpubr::stat_compare_means(comparisons = list(c('II','IT')), vjust = 1.5, method = 't.test')+ expand_limits(y = 0) -> p1

ggplot(total_length, aes(group, IGG, color = group)) +
  geom_jitter(width = 0.05, height = 0, size = 3) +
  stat_summary(geom = 'crossbar', width = 0.3, fun = 'mean', color = 'red') +
  theme_pubr(legend = 'right') +
  labs(title = 'OVA-specific IgG level', y = '') +
  stat_compare_means(comparisons = list(c('II','IT')), vjust = 1.5, method = 't.test')+ expand_limits(y = 0) -> p2

p1 + p2 

# MC38 tumor size plot -------
read_tsv('../CRC-I/data/other/MC38-II-IT-1026.txt') -> mc38

mc38 %>%
  pivot_longer(3:9, names_to = 'day', values_to = 'size') -> mc38

mc38$day <- as.numeric(mc38$day)

mc38 %>%
  ggplot(aes(day, size, color = genotype, group = genotype))+
  stat_summary(fun = 'mean', geom = 'point', size = 2) +
  stat_summary(fun.data = 'mean_se', geom = 'errorbar', width = 1, size = 1) +
  stat_summary(fun = 'mean', geom = 'path', size = 1)+
  ggprism::theme_prism() +
  labs(x = 'days after inoculation', y = 'tumor size (mm3)', title = 'MC38 s.c. tumor growth') +
  stat_compare_means(label = 'p.signif') +
  scale_color_discrete(label = c('GR (n=4)','IT-GR (n=2)'))

