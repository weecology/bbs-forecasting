#!/apps/R/3.2.0/bin/Rscript

devtools::load_all()
library(parallel)
library(tidyverse)
library(gbm)
settings = yaml::yaml.load_file("settings.yaml")

if (!file.exists("observer_model.rds")) {
  fit_observer_model()
} else {
  obs_model = readRDS("observer_model.rds")
}

dir.create("sdm_predictions", showWarnings = FALSE)

settings = yaml::yaml.load_file("settings.yaml")

# GBM parameters
first_validation_year = 2000 # for temporal cross-validation
interaction.depth = 8
n_trees_to_add = 1000
max_trees = 1E4
my_formula = as.formula(paste("present ~", 
                              paste(settings$vars, collapse = "+")))

# Munge -------------------------------------------------------------------
d = get_pop_ts_env_data(settings$start_yr, settings$end_yr, 
                        settings$min_num_yrs) %>% 
  filter(!is.na(species_id))

species_ids = unique(d$species_id)

sp_id = d$species_id[1]
iter = 1

all_data = d %>% 
  filter(species_id == sp_id) %>% 
  select(site_id, year, species_id, observer_id, abundance) %>% 
  right_join(obs_model$data, by = c("site_id", "year", "observer_id")) %>% 
  mutate(present = ifelse(is.na(abundance), 0, 1))

# Fitting/predicting ------------------------------------------------------
fit_species = function(sp_id, iter, all_data){
  cat(as.character(Sys.time()), ": fitting species", sp_id, "\n")
  
  # Training data must be arranged by year for temporal cross-validation
  train = all_data %>% 
    filter(in_train, iteration == iter) %>% 
    arrange(year)
  test = all_data %>% 
    filter(!in_train, iteration == iter)
  
  train_occupancy = train %>% 
    filter(year < first_validation_year) %>% 
    summarize(mean = mean(present), var = var(present))
  
  
  message("think about 0-variance species")
  # avoid problems with 0% or 100% presence
  if (train_occupancy$var == 0) {
    return(
      list(
        predictions = train_occupancy$mean,
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
    data = train,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = mean(train$year < first_validation_year),
    shrinkage = 0.01,
    n.trees = 2 * n_trees_to_add
  )
  
  while (!done_adding) {
    total_trees = probe_model$n.trees
    n.trees = gbm.perf(probe_model, method = "test", plot.it = FALSE)
    
    # We're done adding trees when the optimal number of trees stops going up
    # or when we've reached the maximum number of trees
    if (n.trees < (total_trees - n_trees_to_add) | total_trees >= max_trees) {
      done_adding = TRUE
    } else{
      probe_model = gbm.more(probe_model, n.new.trees = n_trees_to_add)
    }
  }
  
  # Fit the full time series with the same number of trees
  prediction_model = gbm(
    my_formula,
    data = train,
    distribution = "bernoulli",
    interaction.depth = interaction.depth,
    train.fraction = 1,
    shrinkage = 0.01,
    n.trees = n.trees
  )
  
  list(
    train_predictions = predict(probe_model, train, n.trees = n.trees, 
                                type = "response"),
    test_predictions = predict(prediction_model, test, n.trees = n.trees, 
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
