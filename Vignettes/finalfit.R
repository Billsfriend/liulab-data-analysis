library(tidyverse)
library(ggpubr)
library(finalfit)

dependent = "differ.factor"

# Specify explanatory variables of interest
explanatory = c("age", "sex.factor", 
                "extent.factor", "obstruct.factor", 
                "nodes")

view(colon_s)

or_plot()