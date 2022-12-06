# https://rdatatable.gitlab.io/data.table/
# improvement on data.frame!
library(data.table)
#DT = as.data.table(iris)

# FROM[WHERE, SELECT, GROUP BY]
# DT  [i,     j,      by]
# [rows, columns, group by..]

#DT[Petal.Width > 1.0, mean(Petal.Length), by = Species]
#      Species       V1
#1: versicolor 4.362791
#2:  virginica 5.552000
#
# data.table’s functions prefixed with “set” and the operator “:=” 
# work without “<-” to alter data without making copies in 
# memory. E.g., the more efficient “setDT(df)” is analogous to
# “df <- as.data.table(df)”.

#DT[1:5,1:3] # select up-left 5*3 data table

# use some flights data as example
flights <- fread('data/flights14.csv')
dim(flights) # dimension means number of rows (observations) and columns (variables)

DT = data.table(
  ID = c("b","b","b","a","a","c"),
  a = 1:6,
  b = 7:12,
  c = 13:18
) # create a data.table from scratch, one column per arguments

setDT(DT) # directly convert a data frame or list into data table without making copies in RAM

# data.table doesn’t set or use row names, ever.

# subset rows
# Get all the flights with “JFK” as the origin airport in the month of June
ans <- flights[origin == "JFK" & month == 6L]

# Get the first two rows from flights
ans <- flights[1:2]
ans <- flights[,year:day]

# Sort flights first by column origin in ascending order, and then by dest in descending order
ans <- flights[order(origin, -dest)]

# select columns
# Select arr_delay column, but return it as a vector
ans <- flights[, arr_delay]

# Select arr_delay column, but return as a data.table instead
ans <- flights[, list(arr_delay)] # .(arr_delay) is the same as list(arr_delay)

# Select both arr_delay and dep_delay columns
ans <- flights[, .(arr_delay, dep_delay)]

# select both arr_delay and dep_delay columns and rename them to delay_arr and delay_dep
ans <- flights[, .(delay_arr = arr_delay, delay_dep = dep_delay)]

# compute or do in j
# How many trips have had total delay < 0?
ans <- flights[, sum( (arr_delay + dep_delay) < 0 )]

# Calculate the average arrival and departure delay for all flights with “JFK” as the origin airport in the month of June
ans <- flights[origin == 'JFK' & month == 6L, .(m_arr = mean(arr_delay), m_dep = mean(dep_delay))]

# How many trips have been made in 2014 from “JFK” airport in the month of June?
ans <- flights[origin == 'JFK' & month == 6L, length(dest)]
# or using .N which denotes the number of rows in the subset
ans <- flights[origin == "JFK" & month == 6L, .N]

# how can I refer to columns by names in j (like in a data.frame)?
flights[, c('arr_delay','dep_delay')]
# or using .. with a vector contain selected names
# .. means 'up-one-level'; as in unix
select_cols = c("arr_delay", "dep_delay")
flights[ , ..select_cols]
# Select columns named in a variable using with = FALSE
flights[, select_cols, with = 0] # with = FALSE means stop recognizing colnames in expression j as variables

# deselect certain cols
flights[, -c('arr_delay','dep_delay')]
flights[, !c('arr_delay','dep_delay')]
flights[, -(arr_delay:dep_delay)]

# group using by
# How can we get the number of trips corresponding to each origin airport?
flights[, .(.N), by = .(origin)]
flights[, .(.N), by = 'origin']
# When there’s only one column or expression to refer to in j and by, we can drop the .() notation. This is purely for convenience
flights[, .N, by = origin]

# How can we calculate the number of trips for each origin airport for carrier code "AA"?
flights[carrier == 'AA', .N, by = origin]

# How can we get the total number of trips for each origin, dest pair for carrier code "AA"?
flights[carrier == 'AA', .N, by = .(origin, dest)]

# How can we get the average arrival and departure delay for each orig,dest pair for each month for carrier code "AA"?
flights[carrier == 'AA', .(mean(arr_delay), mean(dep_delay)), by = .(origin, dest,month)]

# Sorted by: keyby
# So how can we directly order by all the grouping variables?
flights[carrier == 'AA', .(mean(arr_delay), mean(dep_delay)), keyby = .(origin, dest,month)]
# keyby is faster than by

# chaining
# How can we order ans using the columns origin in ascending order, and dest in descending order?
flights[carrier == "AA", .N, by = .(origin, dest)][order(origin, -dest)]

# Expressions in by
# how many flights started late but arrived early (or on time), started and arrived late etc…
flights[, .N, .(dep_delayed = dep_delay>0, arr_delayed = arr_delay>0)]
# The last row corresponds to dep_delay > 0 = TRUE and arr_delay > 0 = FALSE

# Multiple columns in j - .SD
# Do we have to compute mean() for each column individually?
# What if you had 100 columns to average mean()?
# use lapply()! [list apply]
# 
# Special symbol .SD: stands for Subset of Data. It by itself is a data.table that holds the data for the current group defined using by
DT
DT[, print(.SD), by = ID]
# .SD contains all the columns except the grouping columns by default.
# To compute on (multiple) columns, we can then simply use the base R function lapply().
DT[, lapply(.SD, mean), by = ID]

# How can we specify just the columns we would like to compute the mean() on?
# .SDcols
flights[carrier == "AA",                       ## Only on trips with carrier "AA"
        lapply(.SD, mean),                     ## compute the mean
        by = .(origin, dest, month),           ## for every 'origin,dest,month'
        .SDcols = c("arr_delay", "dep_delay")] ## for just those specified in .SDcols

# Subset .SD for each group:
# How can we return the first two rows for each month?
ans <- flights[, head(.SD, 2), by = month]

# j is very flexiable
# How can we concatenate columns a and b for each group in ID?
DT[, .(val = c(a,b)), by = ID] # basic function c() [concatenate]

# What if we would like to have all the values of column a and b concatenated, but returned as a list column?
DT[, .(val = list(c(a,b))), by = ID]
# A list column can contain any object in each cell, and in this example, each cell is itself a vector and some cells contain longer vectors than others.

# And remember the tip:
# As long as j returns a list, each element of the list will become a column in the resulting data.table.