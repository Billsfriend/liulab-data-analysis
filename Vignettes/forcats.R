library(tidyverse)

f <- factor(c('a','c','b','c'))
f

levels(f)

table(f)
fct_count(f)

fct_unique(f)

# change order of levels
# mentioned level(s) will be put in front
fct_relevel(f, c('c'))

# order levels by freq
fct_infreq(c('a','c','c'))

# order numeric levels by value
fct_inseq(as_factor(c(9,4,-1)))

# order by the order in which they first appear
fct_inorder(c('b_cell','another','t cell'))

# move an end to another end
fct_shift(factor(c('a','b','c')))


ggplot(iris, aes(Species, Sepal.Width, colour = Species)) +
  geom_boxplot()

# fct_reorder() change order by another variable
# Note that lines match order in legend
iris %>%
  mutate(Species = fct_reorder(Species, Sepal.Width)) %>%
  ggplot(aes(Species, Sepal.Width, colour = Species)) +
  geom_boxplot()

chks <- subset(ChickWeight, as.integer(Chick) < 10)
chks <- transform(chks, Chick = fct_shuffle(Chick))

ggplot(chks, aes(Time, weight, colour = Chick)) +
  geom_point() +
  geom_line()

# reorder by final value plot by 2 var
ggplot(chks, aes(Time, weight, colour = fct_reorder2(Chick, Time, weight))) +
  geom_point() +
  geom_line() +
  labs(colour = "Chick")

# change some level(s) to 'other'
fct_other(f, keep = c("a", "b"))
# collapse levels appearing fewer than 2 times into 'other'
fct_lump_min(f, min = 2)
