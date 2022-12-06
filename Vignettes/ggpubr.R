library(ggpubr)
library(tidyverse)

data('ToothGrowth')
df <- ToothGrowth

# boxplot with jitter points ------
ggboxplot(df, x = 'dose', y = 'len', color = 'dose', add = 'jitter', shape = 'dose') ->
  p

# add a layer of p-value
my_comparing <- list(c('0.5','1'), c('0.5','2'))

p + stat_compare_means(comparisons = my_comparing)

# add global p value
p + stat_compare_means(label.y = 50)

# violin plot with box plot
ggviolin(df,
         x = "dose", y = "len", fill = "dose",
         palette = c("#00AFBB", "#E7B800", "#FC4E07"),
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparing, label = "p.signif") # Add significance levels ***
