library(tidyverse)

browseVignettes('dplyr') # open vignettes in browser

dim(starwars)
# a tibble: a improved class of data.frame

# filter() select observation by variables
starwars %>% filter(skin_color == "light", eye_color == "brown")
# use pipe operator %>% 'then'

# equivalent to this base R code
starwars[starwars$skin_color == "light" & starwars$eye_color == "brown", ]

# arrange() reorder observation by variables
starwars %>% arrange(height, mass)

# desc() means descending order
starwars %>% arrange(desc(height))
