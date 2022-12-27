# 2022/4/29

library(data.table)
library(tidyverse)

data <- fread('data/sumpercent.txt')

# barplot with errorbar
ggplot(data, aes(x=subset, y=mean, fill=genotype))+
  geom_col(position = 'dodge')+
  geom_errorbar(aes(ymin = mean - CI, ymax = mean + CI), position = position_dodge(0.9), width = .2)+
  ylab('Proportion in all blood immune cells')+
  theme_prism(base_size = 14)


# t-test
iisubset <- 552
iiall <- 1468
itsubset <- 3275
itall <- 9226

ii <- c(rep(1, iisubset), rep(0, iiall-iisubset))
it <- c(rep(1, itsubset), rep(0, itall-itsubset))
# two sample t-test
t.test(ii, it)
