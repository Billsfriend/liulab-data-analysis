library(dplyr)

starwars %>%
  summarise(across(where(is.character), n_distinct))

starwars %>%
  group_by(species) %>%
  filter(n() > 1) %>%
  summarise(across(c(sex, gender, homeworld), n_distinct))

starwars %>% 
  group_by(homeworld) %>% 
  filter(n() > 1) %>% 
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))

# doesn’t select grouping variables in order to avoid accidentally modifying them
df <- data.frame(g = c(1, 1, 2), x = c(-1, 1, 3), y = c(-1, -4, -9))
df %>% 
  group_by(g) %>% 
  summarise(across(where(is.numeric), sum))

# multiple functions ---------
# a named list of functions
min_max <- list(
  min = ~min(.x, na.rm = TRUE),
  max = ~max(.x, na.rm = TRUE)
)

starwars %>%
  summarise(across(where(is.numeric), min_max))

# define names of summarized columns by .name and glue spec
starwars %>%
  summarise(across(where(is.numeric), min_max, .names = '{.fn}.{.col}'))

# current column -----------
# get name of current column by cur_column()
df <- tibble(x = 1:3, y = 3:5, z = 5:7)
mult <- list(x = 1, y = 10, z = 100)

df %>%
  mutate(across(all_of(names(mult)), ~ .x * mult[[cur_column()]]))

# gotchas --------
df <- data.frame(x = c(1, 2, 3), y = c(1, 4, 9))

# summarise() compute each variable step-wise
df %>%
  summarise(n = n(), across(where(is.numeric), sd))

df %>%
  summarise(across(where(is.numeric), sd), n = n())

df %>% 
  summarise(n = n(), across(where(is.numeric) & !n, sd))

# or define a new tibble so the result variable do not affect next step input
df %>% 
  summarise(
    tibble(n = n(), across(where(is.numeric), sd))
  )

# other verbs ---------
rescale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

df <- tibble(x = 1:4, y = rnorm(4))
df %>%
  mutate(across(where(is.numeric), rescale01))

# for some verbs, the .fun can be omitted
starwars %>% distinct(across(contains('color')))

starwars %>% count(across(contains('color'), sort = TRUE))

# filter() ---------
# if_any() keep the rows for at least one TRUE column
starwars %>%
  filter(if_any(everything(), ~ !is.na(.x)))

# if_all() keep the row if all col is TRUE
starwars %>%
  filter(if_all(everything(), ~ !is.na(.x)))

# this will throw a warning
starwars %>%
  filter(across(everything(), ~ !is.na(.x)))


# ROWWISE ----------
# rowwise() is kind of a special group_by() that every row is a group
df <- tibble(x = 1:2, y = 3:4, z = 5:6)

df %>% mutate(m = mean(c(x, y, z)))
df %>% rowwise() %>% mutate(m = mean(c(x, y, z)))

# pass a row identifier
df <- tibble(name = c("Mara", "Hadley"), x = 1:2, y = 3:4, z = 5:6)
df %>% 
  rowwise(name) %>% 
  summarise(m = mean(c(x, y, z)))


df <- tibble(id = 1:6, w = 10:15, x = 20:25, y = 30:35, z = 40:45)

# use c_across() as c() with tidy-select
rf <- rowwise(df)

rf %>% mutate(total = sum(c_across(where(is.numeric))))

rf %>% 
  mutate(total = sum(c_across(w:z))) %>% 
  ungroup() %>% 
  mutate(across(w:z, ~ . / total))

# row-wise summary functions ---------
df %>% mutate(total = rowSums(across(where(is.numeric))))

# list-columns --------
df <- tibble(x = list(1, 2:3, 4:6))

df %>% mutate(l = sapply(x, length))
df %>% rowwise() %>% mutate(l = length(x))

# key difference is that when mutate() slices up the columns to pass to length(y) the grouped mutate uses [ and the row-wise mutate uses [[.

# modeling ---------
# nest_by() create nested tibble that is auto-row-wise
by_cyl <- mtcars %>% nest_by(cyl)
by_cyl

mods <- by_cyl %>% mutate(mod = list(lm(mpg ~ wt, data = data)))
mods

mods <- mods %>% mutate(pred = list(predict(mod, data)))

mods %>% summarise(rmse = sqrt(mean((pred - data$mpg) ^ 2)))

mods %>% summarise(broom::tidy(mod))

# repeated function calls ----------
df <- tribble(
  ~n, ~min, ~max,
  1, 0, 1,
  2, 10, 100,
  3, 100, 1000,
)

# generate random deviates from uniform distribution
df %>%
  rowwise() %>%
  mutate(data = list(runif(n, min, max)))

# multiple combinations ------------
# expand.grid() generates combination in data frame
df <- expand.grid(mean = c(-1, 0, 1), sd = c(1, 10, 100))

df %>%
  rowwise() %>%
  mutate(data = list(rnorm(10, mean, sd)))

# varying functions -------
df <- tribble(
  ~rng, ~params,
  'runif', list(n = 10),
  'rnorm', list(n = 20),
  'rpois', list(n = 10, lambda = 5),
) %>%
  rowwise()

# do.call() invoke by function names and params
df %>%
  mutate(data = list(do.call(rng, params)))


