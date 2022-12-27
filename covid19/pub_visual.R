library(tidyverse)
library(ggpubr)

read_tsv('covid19/data/odds_ratio_forest_plot.tsv') -> frst_data

# Korea cohort --------
frst_data %>%
  filter(str_detect(cohort, 'Korea')) %>%
  ggplot(aes(model, y = `Odds Ratio`))+
  geom_pointrange(shape = 'square', color = 'blue', size = 2,aes(ymax = L95, ymin = U95))+
  geom_text(aes(label = `Odds Ratio`), nudge_x = 0.4)+
  scale_y_log10(breaks = c(0.4,0.7,1,1.3,1.6))+
  theme_pubr()+
  labs_pubr()+
  coord_flip()+
  facet_wrap(~comparison, ncol = 1)+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  expand_limits(y = c(0.3,1.7)) +
  ggtitle('South Korea cohort (6500 HC, 108 cases, 69 hospitalized)')

last_plot() -> forest_plot_korea

frst_data %>%
  mutate(`Odds Ratio` = as.character(`Odds Ratio`)) %>%
  mutate(`P-value` = as.character(`P-value`)) %>%
  pivot_longer(cols = c(`Odds Ratio`, `P-value`, `95% CI`), names_to = 'info', values_to = 'value') %>%
  mutate(info = fct_relevel(as_factor(info), c('Odds Ratio','95% CI','P-value'))) %>%
  filter(str_detect(cohort, 'Kore')) %>%
  ggplot(aes(y = model, x = info, label = value))+
  scale_x_discrete(position = 'top')+
  geom_text()+
  theme_classic()+
  facet_wrap(~comparison, ncol = 1)+
  theme(
    axis.line = element_blank(),
    axis.ticks  = element_blank(),
    axis.title.y  = element_blank(),
    axis.title.x  = element_blank(),
    axis.text.x = element_text(color="white")
  )+
  labs_pubr()

last_plot() -> data_table_korea

forest_plot_korea + data_table_korea

# Japan cohort
frst_data %>%
  filter(str_detect(cohort, 'Jap')) %>%
  mutate(`P-value` = str_sub(str_glue('p={`P-value`}'), end = 7)) %>%
  ggplot(aes(model, y = `Odds Ratio`))+
  geom_pointrange(shape = 'square', color = 'blue', size = 2,aes(ymax = L95, ymin = U95))+
  geom_text(aes(label = `Odds Ratio`), nudge_x = 0.3)+
  geom_text(aes(label = `P-value`, y = 1.4))+
  scale_y_log10(breaks = c(0.4,0.7,1,1.3,1.6))+
  theme_pubr()+
  labs_pubr()+
  coord_flip()+
  facet_wrap(~comparison, scales = 'free_x', ncol = 1)+
  geom_hline(yintercept = 1, linetype = 'dashed')+
  expand_limits(y = c(0.6,1.7)) +
  ggtitle('Japan cohort (3289 HC, 2393 cases, 1980 severe cases)')+
  xlab('grouping')

last_plot() -> p2

p1 / p2

ggsave('covid19/figures/EastAsia_forest_plot.pdf',
       width = 24,
       height = 16,
       units = 'cm')

# Wuhan cohort ---------
frst_data %>%
  filter(str_detect(cohort, 'Wuhan')) %>%
  mutate(`P-value` = str_sub(str_glue('p={`P-value`}'), end = 7)) %>%
  ggplot(aes(comparison, y = `Odds Ratio`))+
  geom_pointrange(shape = 'square', color = 'blue', size = 2,aes(ymax = L95, ymin = U95))+
  geom_text(aes(label = `Odds Ratio`), nudge_x = 0.3)+
  geom_text(aes(label = `P-value`, y = 1.8), nudge_x = 0.2)+
  scale_y_log10()+
  theme_pubr()+
  labs_pubr()+
  coord_flip()+
  geom_hline(yintercept = 1, linetype = 'dashed') +
  ggtitle('Wuhan cohort (163 cases)')+
  xlab('grouping')

frst_data %>%
  filter(str_detect(cohort, 'Shenzhen')) %>%
  mutate(`P-value` = str_sub(str_glue('p={`P-value`}'), end = 7)) %>%
  mutate(comparison = fct_relevel(comparison, 'Additive','Dominative', 'Recessive')) %>%
  ggplot(aes(comparison, y = `Odds Ratio`))+
  geom_pointrange(shape = 'square', color = 'blue', size = 2,aes(ymax = L95, ymin = U95))+
  geom_text(aes(label = str_sub(`Odds Ratio`, end = 5)), nudge_x = 0.3)+
  geom_text(aes(label = `P-value`, y = 5), nudge_x = 0.2)+
  scale_y_log10()+
  theme_pubr()+
  labs_pubr()+
  coord_flip()+
  geom_hline(yintercept = 1, linetype = 'dashed', alpha = 0.3) +
  ggtitle('Shenzhen cohort (233 Cases, 459 HC)')+
  xlab('Model')

# HGI data -----------
hgi_my_snp <- read_csv('covid19/results/hgi_FCGR2B_I232T.csv')

hgi_my_snp %>%
  mutate(odds_ratio = exp(all_inv_var_meta_beta),
         U95 = exp(all_inv_var_meta_beta + 1.96*all_inv_var_meta_sebeta),
         L95 = exp(all_inv_var_meta_beta - 1.96*all_inv_var_meta_sebeta),
         p_value = str_sub(str_glue('p={all_inv_var_meta_p}'), end = 7),
         population = case_when(
           popu == 'afr' ~ 'Africa',
           popu == 'all' ~ 'Global',
           popu == 'eas' ~ 'East Asia',
           popu == 'eur' ~ 'Europe'
         ),
         group = case_when(
           group == 'B2' ~ 'Hospitalized cases vs HC',
           group == 'C2' ~ 'All cases vs HC',
           TRUE ~ group
         ),
         size = str_glue('{all_inv_var_meta_cases} cases vs {all_inv_var_meta_controls} HC')
         ) ->
  hgi_my_snp

write_csv(hgi_my_snp, 'covid19/results/hgi_FCGR2B_I232T.csv')

hgi_my_snp %>%
  ggplot(aes(population, y = odds_ratio))+
  geom_pointrange(shape = 'square', color = 'blue', size = 2,aes(ymax = L95, ymin = U95))+
  geom_text(aes(label = str_sub(odds_ratio, end = 5)), nudge_x = 0.3)+
  geom_text(aes(label = p_value, y = 1.3), nudge_x = 0.2)+
  geom_text(aes(label = size, y = 1.3), nudge_x = -0.2)+
  scale_y_log10()+
  theme_pubr()+
  labs_pubr()+
  coord_flip()+
  facet_wrap(~group)+
  geom_hline(yintercept = 1, linetype = 'dashed', alpha = 0.3) +
  expand_limits(y = 1.5)+
  ggtitle('HGI data')+
  xlab('Cohort')
