devtools::load_all()
library(tidyverse)
library(rstan)
start_yr = 1982
end_yr = 2013
min_num_yrs = 25
last_train_year = 2003
d = get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>% 
  filter(!is.na(observer_id))
train = d[d$year <= last_train_year, ]
test = d[d$year > last_train_year, ]

data = list(
  N = nrow(train),
  N_observers = length(unique(train$observer_id)),
  N_sites = length(unique(train$site_id)),
  y = c(scale(train$richness)),
  site_index = as.integer(factor(train$site_id)),
  observer_index = as.integer(factor(train$observer_id))
)

s = stan("lme.stan", data = data, cores = 4)

saveRDS(s, file = "lme-stan.rds")
