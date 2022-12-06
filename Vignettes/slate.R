library(tidyverse)
library(ggprism)
# browseVignettes('ggprism')

data <- read_tsv('data/other/220625-qP-RelaExp.txt')

data %>%
  filter(Time == '2') %>%
  ggplot(aes(x = Substance, y = RelativeExp)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'IL-1B 2h')

data %>%
  filter(Time == '6') %>%
  ggplot(aes(x = Substance, y = RelativeExp)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  coord_flip()+
  labs(title = 'IL-1B 6h')

data %>%
  filter(Gene == 'IL-10' & Group == '6h') %>%
  ggplot(aes(x = Name, y = fc)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'IL-10 6h')

data %>%
  filter(Gene == 'IL-1B' & Group == '2h') %>%
  ggplot(aes(x = Name, y = fc)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'IL-1B 2h')

data %>%
  filter(Gene == 'IL-1B' & Group == '6h') %>%
  ggplot(aes(x = Name, y = fc)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'IL-1B 6h')

data %>%
  filter(Gene == 'IL-1B') %>%
  ggplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', aes(x = Name, y = fc, fill = Group)) +
  labs(title = 'IL-1B')

data %>%
  filter(Gene == 'TNFA' & Group == '2h') %>%
  ggplot(aes(x = Name, y = fc)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'TNFA 2h')

data %>%
  filter(Gene == 'TNFA' & Group == '6h') %>%
  ggplot(aes(x = Name, y = fc)) +
  geom_dotplot(binaxis = 'y', stackdir = 'center') +
  stat_summary(geom = "crossbar", fun = mean, colour = "red") +
  ylim(0,8)+
  labs(title = 'TNFA 6h')

e <- ggplot(mpg, aes(x=cty, y=hwy))
e + geom_point()
e + geom_quantile()
e + geom_smooth()
e + stat_ecdf(n=40)
e + stat_boxplot(coef = 1.5)
e + stat_sum()
e + stat_summary(fun.data = 'mean_cl_boot')
e + stat_identity()
e + stat_unique()

data %>%
  filter(Gene == 'TNFA' & Group == '6h') %>%
  ggplot() +
  geom_boxplot(aes(x = Name, y = fc)) +
  labs(title = 'TNFA 6h')
