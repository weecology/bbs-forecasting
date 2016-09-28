#!/apps/R/3.2.0/bin/Rscript

devtools::load_all()
library(parallel)
library(tidyverse)
library(gbm)

dir.create("sdm_predictions", showWarnings = FALSE)

start_yr = 1982
end_yr = 2013
min_num_yrs = 25
last_train_year = 2003

# GBM parameters
first_validation_year = 2000 # for temporal cross-validation
interaction.depth = 8
max_n_trees = 2E4
my_formula = present ~ bio2 + bio5 + bio15 + ndvi_sum + elevs +
  observer_effect + site_effect

# Munge -------------------------------------------------------------------
d = get_pop_ts_env_data(start_yr, end_yr, min_num_yrs) %>% 
  add_ranefs(last_training_year = last_train_year)

train = d[d$year <= last_train_year, ]
test = d[d$year > last_train_year, ]

# Training data must be arranged by year first for temporal cross-validation
# to work!
train_x = select(train, -species_id, -abundance) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  distinct() %>% 
  arrange(year, site_id)

test_x = select(test, -species_id, -abundance) %>% 
  distinct() %>% 
  filter(!is.na(lat)) %>% 
  distinct() %>% 
  arrange(year, site_id)

species_ids = unique(train$species_id)


# Fitting/predicting ------------------------------------------------------
fit_species = function(sp_id){
  cat(strsplit(date(), " ")[[1]][4], " fitting species", sp_id, "\n")
  
  xy = filter(train, species_id == sp_id) %>% 
    select(site_id, year) %>% 
    mutate(present = 1) %>% 
    right_join(train_x, c("site_id", "year")) %>% 
    arrange(year, site_id) %>% 
    mutate(present = ifelse(is.na(present), 0, 1))
  
  # avoid problems with 100% or 0% presence
  if (var(xy$present[xy$year < first_validation_year]) == 0) {
    return(NA)
  }
  
  
  # Determine the number of trees to include by using the first portion of
  # the training set to predict the remainder of the training set.
  g = gbm(
    my_formula,
    data = xy,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = mean(xy$year < first_validation_year),
    n.trees = max_n_trees
  )
  n.trees = gbm.perf(g)
  
  # Fit the full time series with the same number of trees
  g_full = gbm(
    my_formula,
    data = xy,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = 1,
    n.trees = n.trees
  )
  
  p = predict(g_full, test_x, n.trees = n.trees, type = "response")
  
  list(
    predictions = p,
    importances = summary(g_full, plot = FALSE),
    n.trees = n.trees
  )
}

results = mclapply(species_ids, fit_species, 
                       mc.cores = 8,
                       mc.preschedule = FALSE)

saveRDS(results, file = "sdm_predictions/gbm_species_results.rds")
