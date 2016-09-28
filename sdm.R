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
n_trees_to_add = 1000
max_trees = 2E4
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
  
  # avoid problems with 0% or 100% presence
  if (mean(xy$present[xy$year < first_validation_year]) == 0) {
    return(
      list(
        predictions = 0,
        importances = NA,
        n.trees = NA
      )
    )
  }
  if (mean(xy$present[xy$year < first_validation_year]) == 1) {
    return(
      list(
        predictions = 1,
        importances = NA,
        n.trees = NA
      )
    )
  }
  
  # Determine the number of trees to include by using the first portion of
  # the training set to predict the remainder of the training set.
  done_adding = FALSE # flag for when we have enough trees
  probe_model = gbm(
    my_formula,
    data = xy,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = mean(xy$year < first_validation_year),
    n.trees = 2 * n_trees_to_add
  )
  
  while (!done_adding) {
    n.trees = gbm.perf(probe_model, method = "test", plot.it = FALSE)
    
    # We're done adding trees when the optimal number of trees stops going up.
    if (n.trees < (probe_model$n.trees - n_trees_to_add) | n.trees > max_trees) {
      done_adding = TRUE
    } else{
      probe_model = gbm.more(probe_model, n.new.trees = n_trees_to_add)
    }
  }
  
  # Fit the full time series with the same number of trees
  prediction_model = gbm(
    my_formula,
    data = xy,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = 1,
    n.trees = n.trees
  )
  
  list(
    train_predictions = predict(probe_model, xy, n.trees = n.trees, 
                                type = "response"),
    test_predictions = predict(prediction_model, test_x, n.trees = n.trees, 
                               type = "response"),
    importances = summary(prediction_model, plot = FALSE),
    n.trees = n.trees,
    valid.error = probe_model$valid.error
  )
}

results = mclapply(species_ids, fit_species, 
                       mc.cores = 8,
                       mc.preschedule = FALSE)

saveRDS(results, file = "sdm_predictions/gbm_species_results.rds")
