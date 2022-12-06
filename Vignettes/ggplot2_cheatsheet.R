library(tidyverse)

# graphical primitives -----
a <- ggplot(economics, aes(date, unemploy))

# show only data range
a + geom_blank()

# x - y as a single path
a + geom_path(lineend = 'square', linejoin = 'round', linemitre = 10)

# connect (x,y) to a polygon
a + geom_polygon(aes(alpha = 50))

# can be used to plot confidence interval
a + geom_ribbon(aes(ymin = unemploy - 900, ymax = unemploy + 900))

b <- ggplot(seals, aes(x = long, y = lat))

b + geom_curve(aes(yend = lat + 1, xend = long + 1), curvature = -0.1)

# draw a rectangle for each row
b + geom_rect(aes(xmin = long, ymin = lat, xmax = long + 0.2, ymax = lat + 0.1))

# plot y = kx + b
b + geom_abline(aes(intercept = 0, slope = 1), color = 'red')

# horizontal and vertical line, mostly as annotation
b + geom_hline(aes(yintercept = lat), linetype = 'dashed')
b + geom_vline(aes(xintercept = long), linetype = 'dotted')

# plot a line with definite begin and end
b + geom_segment(aes(yend = lat + 1/2, xend = long + 1/2))

# a line with a angle as direction and radius as length
b + geom_spoke(aes(angle = 1:1155, radius = 1))

# one variable continuous -------
c <- ggplot(mpg, aes(hwy)) ; c2 <- ggplot(mpg)

# distribution with auto-selected bin
c + geom_area(stat = 'bin')

# smoothed distribution function
c + geom_density()

# dot plot as distribution
c + geom_dotplot()

# like geom_area with no fill
c + geom_freqpoly()

# barplot-like histogram
c + geom_histogram()

# quantile-quantile plot to test if sample fit with a certain distribution
c2 + geom_qq(aes(sample = hwy))

# one var discrete -----
d <- ggplot(mpg, aes(fl))

# count each discrete elements
d + geom_bar()

# two var both continuous ------
e <- ggplot(mpg, aes(cty, hwy))

# add label with text for each row
e + geom_label(aes(label = cty), nudge_x = 1, nudge_y = 1)

# add text without label box for each row
e + geom_text(aes(label = cty), nudge_x = 1, nudge_y = 1)

# classic x,y point
e + geom_point()

# fits a quantile regression and draw fit lines
e + geom_quantile()

# rug plot: compact visualisation designed to supplement a 2d display with the two 1d marginal distributions
# best used with smaller datasets
e + geom_rug()

# smooth fit formula
e + geom_smooth(method = lm)

# randomly moved geom_point
e + geom_jitter(height = 2, width = 2)

# one discrete, one continuous --------
f <- ggplot(mpg, aes(class, hwy))

# column plot with numeric data
f + geom_col()

# boxplot with default quantile
f + geom_boxplot()

# dot plot like one var distribution, but on y-axis of each x var
f + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3)

# violin plot.
# scale = 'area'(default) means every violin has the same area. 'count' will reflect count number of each group on area
f + geom_violin(scale = 'area')

# both discrete ------
g <- ggplot(diamonds, aes(cut, color))

# count for every x-y combination
g + geom_count()

# continuous bivariate distribution -----
h <- ggplot(diamonds, aes(carat, price))

# 2d density distribution in matrix
h + geom_bin2d(binwidth = c(0.25, 500))

# contour line density plot
h + geom_density_2d()

# hexagon density
h + geom_hex()

# continuous function ------
a + geom_area()

# much like geom_path?
a + geom_line()

# make line steppy
a + geom_step(direction = 'hv')

# visualize error -----
df <- data.frame(grp = c('A','B'), fit = 4:5, se = 1:2)

j <- ggplot(df, aes(grp, fit, ymin = fit - se, ymax = fit + se))

# box-like error bar
j + geom_crossbar()

# T-shape error bar
j + geom_errorbar()

# just vertical line...
j + geom_linerange()

# big point at center
j + geom_pointrange()

# maps ------
data <- data.frame(murder = USArrests$Murder, state = tolower(rownames(USArrests)))

map <- map_data('state')

k <- ggplot(data, aes(fill = murder))

k + geom_map(aes(map_id = state), map = map) +
  expand_limits(x = map$long, y = map$lat)

# three variables -------
seals$z <- with(seals, sqrt(delta_long^2 + delta_lat^2))

l <- ggplot(seals, aes(x = long, y = lat, z = z))

# visualize 3d height (z) in contour
l + geom_contour()

# fill with color
l + geom_contour_filled()

# matrix like
l + geom_tile(aes(fill = z))

# geom_raster() is a high performance special case for when all the tiles are the same size
l + geom_raster(aes(fill = z), hjust = 0.5, vjust = 0.5, interpolate = FALSE)

# stats layer --------
a <- ggplot(economics, aes(date, unemploy))

a + stat_density_2d(aes(fill = ..level..), geom = 'polygon')

c <- ggplot(mpg, aes(hwy))

# one var
c + stat_bin(binwidth = 1, boundary = 10)

c + stat_count(width = 1)

c + stat_density(adjust = 1, kernel = 'gaussian')

# two conti var
e <- ggplot(mpg, aes(cty, hwy))

e + stat_bin_2d(bins = 30, drop = TRUE)

e + stat_bin_hex(bins = 30)

e + stat_density_2d(contour = TRUE, n = 100)

# compute normal data ellipse
e + stat_ellipse(level = 0.95, segments = 51, type = 't')

l <- ggplot(seals, aes(long, lat))

l + stat_contour(aes(z = z))

l + stat_summary_hex(aes(z = z),bins = 30, fun = max)

l + stat_summary_2d(aes(z = z),bins = 30, fun = mean)

f <- ggplot(mpg, aes(class, hwy))

f + stat_boxplot(coef = 1.5)

# like violin plot
f + stat_ydensity(kernel = 'gaussian', scale = 'area')

e <- ggplot(mpg, aes(cty, hwy))

e + stat_ecdf(n = 40)

e + stat_quantile(quantiles = c(0.1, 0.9), formula = y ~ log(x), method = 'rq')

e + stat_smooth(method = 'lm', formula = y ~ x, se = TRUE, level = 0.95)

e + stat_sum()

e + stat_summary(fun.data = 'mean_cl_boot')

e + stat_identity()

e + stat_unique()

h <- ggplot(diamonds, aes(carat, price))

h + stat_summary_bin(fun = 'mean', geom = 'bar')

# plot a function
ggplot() + xlim(-5, 5) + stat_function(fun = dnorm, n = 20, geom = 'point')

ggplot() + stat_qq(aes(sample = 1:100))

# scales package -----
d <- ggplot(mpg, aes(fl))
n <- d + geom_bar(aes(fill = fl))

n + scale_fill_manual(values = c('skyblue', 'royalblue', 'blue', 'navy'),
                      limits = c('d', 'e', 'p', 'r'),
                      breaks = c('d', 'e', 'p', 'r'),
                      name = 'fuel',
                      labels = c('D', 'E', 'P', 'R'))

# log axis
n + scale_y_log10()

# reverse axis
n + scale_y_reverse()

# square root scale
n + scale_y_sqrt()

# discrete color
RColorBrewer::display.brewer.all()

n + scale_fill_brewer(palette = 'Blues')

n + scale_fill_grey(start = 0.2, end = 0.8, na.value = 'red')

# disable legends
n + guides(fill = 'none')

# change position of legends
n + theme(legend.position = "bottom")

# set legend title and labels
n + scale_fill_discrete(name = "Title", 
                        labels = c("A", "B", "C", "D", "E")) 

# continuous color
o <- c + geom_dotplot(aes(fill = ..x..))

o + scale_fill_distiller(palette = 'Blues')

# specify two color into gradient
o + scale_fill_gradient(low = 'red', high = 'yellow')

# specify three color and midpoint
o + scale_fill_gradient2(low = 'red', high = 'blue', mid = 'white', midpoint = 25)

# n color as gradient
o + scale_fill_gradientn(colors = topo.colors(6))

# palettes
rainbow(7)
heat.colors(9)
terrain.colors(5)
cm.colors(6)

# shape and size ----
e <- ggplot(mpg, aes(cty, hwy))

p <- e + geom_point(aes(shape = fl, size = cyl))

p + scale_shape() + scale_size()

# check document to list shape with index
p + scale_shape_manual(values = c(3:7))

p + scale_radius(range = c(1, 6))

p + scale_size_area(max_size = 6)

# coordinate systems
r <- d + geom_bar()

r + coord_cartesian(xlim = c(0, 5))

r + coord_fixed(ratio = 1/3)

r + coord_polar(theta = 'x', direction = 1)

r + coord_trans(y = 'sqrt')

# map coordinate
worldmap <- map_data('world')

ggmap <- ggplot(worldmap, aes(x = long, y = lat, group = group)) +
  geom_polygon(fill = "white", colour = "black")

ggmap

ggmap + coord_quickmap()

# facet ---------
t <- ggplot(mpg, aes(cty, hwy)) + geom_point()

t + facet_grid(cols = vars(fl))
t + facet_grid(cols = vars(fl), labeller = label_both)

t + facet_grid(rows = vars(year))

t + facet_grid(cols = vars(fl), rows = vars(year))

# wrap into rectangle layout
t + facet_wrap(vars(fl))

# set scales to let axis vary across facets
t + facet_grid(cols = vars(fl), rows = vars(drv), scales = 'free')

# labels and legends --------
t + labs(x = 'x axis label',
         y = 'y axis label',
         title = 'title above plot',
         subtitle = 'subtitle below title',
         caption = 'caption below plot')

# manually place a geom
t + annotate(geom = "text", x = 8, y = 9, label = 'A')

