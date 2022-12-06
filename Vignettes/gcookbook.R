# https://r-graphics.org/
# install.packages('gcookbook')

library(tidyverse)
library(gcookbook)

simpledat
#>    A1 A2 A3
#> B1 10  7 12
#> B2  9 11  6

# barplot plot a column as a group as default
barplot(simpledat, beside = TRUE)

# transpose
barplot(t(simpledat), beside=TRUE)

# plot line, separately
plot(simpledat[1,], type="l")
lines(simpledat[2,], type="l", col="blue")

# With ggplot2, the structure of the data is always the same: it requires a data frame in “long” format, as opposed to the “wide” format used previously
# here is long format of the data used previously
simpledat_long

ggplot(simpledat_long, aes(x = Aval, y = value, fill = Bval)) +
  geom_col(position = "dodge")

# no need to transpose data, just exchange x and fill
ggplot(simpledat_long, aes(x = Bval, y = value, fill = Aval)) +
  geom_col(position = "dodge")

# change geom_col to geom_line. and fill to color
ggplot(simpledat_long, aes(x = Aval, y = value, colour = Bval, group = Bval)) +
  geom_line()

# make up a data frame
dat <- data.frame(
  xval = 1:4,
  yval=c(3, 5, 6, 9),
  group=c("A","B","A","B")
)

dat
#>   xval yval group
#> 1    1    3     A
#> 2    2    5     B
#> 3    3    6     A
#> 4    4    9     B

# plot point
ggplot(dat, aes(x = xval, y = yval)) +geom_point()

# plot specification can be stored and reused
p <- ggplot(dat, aes(x = xval, y = yval))

p + geom_point()
p + geom_point(aes(color = group)) # add color mapping by group
p + geom_point(color = 'blue') # not aesthetic mapping, just set color

p +
  geom_point() +
  scale_x_continuous(limits = c(0, 8)) # change the scale of x

p +
  geom_point(aes(colour = group)) +
  scale_colour_manual(values = c("orange", "forestgreen")) # modify the default color scale

# the pipe operator, %>%
library(dplyr) # The pipe is provided by dplyr

morley # Look at the morley data set
#>     Expt Run Speed
#> 001    1   1   850
#> 002    1   2   740
#> 003    1   3   900
#>  ...<94 more rows>...
#> 098    5  18   800
#> 099    5  19   810
#> 100    5  20   870

morley %>%
  filter(Expt == 1) %>%
  summary()
#>       Expt        Run            Speed     
#>  Min.   :1   Min.   : 1.00   Min.   : 650  
#>  1st Qu.:1   1st Qu.: 5.75   1st Qu.: 850  
#>  Median :1   Median :10.50   Median : 940  
#>  Mean   :1   Mean   :10.50   Mean   : 909  
#>  3rd Qu.:1   3rd Qu.:15.25   3rd Qu.: 980  
#>  Max.   :1   Max.   :20.00   Max.   :1070

# without the pipe operator, the code would like this:
summary(filter(morley, Expt == 1))

# %>% connect the function and make them more readable

# h(g(f(x)))
# 
# # Equivalent to:
# x %>%
#   f() %>%
#   g() %>%
#   h()

# If you want to store the final result, you can use the <- operator at the beginning
# x <- x %>%
#   f() %>%
#   g() %>%
#   h()

# if there is additional argument for function, the pipe will put the value to the right. these 2 line are equivalent:
filter(morley, Expt == 1)

morley %>% filter(Expt == 1)

# boxplot: auto-calculate mean and sd
ggplot(ToothGrowth, aes(x = supp, y = len)) +
  geom_boxplot()
# interaction() to enclose multiple groups
ggplot(ToothGrowth, aes(x = interaction(supp, dose), y = len)) +
  geom_boxplot()

# reorder the data in graph
# Make a copy of the InsectSprays data set since we're modifying it
iss <- InsectSprays
iss$spray
#>  [1] A A A A A A A A A A A A B B B B B B B B B B B B C C C C C C C C C C C C
#> [37] D D D D D D D D D D D D E E E E E E E E E E E E F F F F F F F F F F F F
#> Levels: A B C D E F

iss$spray <- reorder(iss$spray, iss$count, FUN = mean)
iss$spray
#>  [1] A A A A A A A A A A A A B B B B B B B B B B B B C C C C C C C C C C C C
#> [37] D D D D D D D D D D D D E E E E E E E E E E E E F F F F F F F F F F F F
#> attr(,"scores")
#>         A         B         C         D         E         F 
#> 14.500000 15.333333  2.083333  4.916667  3.500000 16.666667 
#> Levels: C E D A B F