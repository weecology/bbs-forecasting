# From forecasting-bbs-R.ipynb --------------------------------------------
library(dplyr)
library(broom)
library(tidyr)
source("forecast-bbs-core.R")

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25


# start SDM ---------------------------------------------------------------
library(gbm)
library(parallel)

mc.cores = 8
n_test_years = 5
interaction.depth = 8
n.trees = 1E4

train_years = start_yr:(end_yr - n_test_years)
test_years = (end_yr - n_test_years + 1):end_yr

# Parameters for spatial cross-validation
n_clusters = 100
n_folds = 10
stopifnot(n_clusters %% n_folds == 0)

# Load data ---------------------------------------------------------------
occurrences = get_bbs_data()
env = get_env_data() %>%
  filter_ts(start_yr, end_yr, min_num_yrs) %>%
  inner_join(distinct(occurrences[ , c("site_id", "lat", "long")]),
             by = "site_id")

# Determine presence-absence ----------------------------------------------

sampling_events = tidyr::expand(occurrences, c(site_id, year)) %>%
  inner_join(env, by = c("site_id", "year"))


make_species_data = function(x){
  out = occurrences %>%
    filter(species_id == x) %>%
    dplyr::select(site_id, year, abundance) %>%
    right_join(sampling_events, by = c("site_id", "year")) %>%
    replace_na(list(abundance = 0))

  if (sum(out$abundance) == 0){
    stop(paste("species", x, "doesn't exist"))
  } else {
    out
  }
}




# Spatial cross validation ------------------------------------------------

set.seed(0)

locations = sampling_events %>%
  distinct(site_id, .keep_all = TRUE) %>%
  dplyr::select(site_id, lat, long)

# assign many clusters of routes that are close in lat/long
cluster_id = kmeans(locations[ , 2:3], n_clusters, nstart = 100)$cluster

# Randomly combine the clusters to make the cross-validation folds
fold_matrix = matrix(sample(1:n_clusters), ncol = n_folds)
fold_id = numeric(nrow(locations))
for (i in 1:n_folds) {
  fold_id = fold_id + i * (cluster_id %in% fold_matrix[, i])
}

locations = cbind(locations, fold_id = fold_id)

# Optionally, view the cross-validation folds
# using `plot(locations[,3:2], col = fold_id, pch = fold_id)`

# fit gbm -----------------------------------------------------------------

# Function to cross-validate a GBM model with different numbers of trees, then
# fit the optimal model size to the full training set.
fit_species = function(species_id){

  data = make_species_data(x = species_id) %>%
    mutate(presence = abundance > 0)

  stopifnot(all(data$lat == sampling_events$lat))

  train_data = data %>%
    filter(year %in% train_years) %>%
    inner_join(locations[ , c("fold_id", "site_id")], by = "site_id") %>%
    dplyr::select(-lat, -long, -site_id, -year, -abundance)

  train_folds = train_data$fold_id

  # Don't want to use fold_id as a predictor variable below
  train_data = dplyr::select(train_data, -fold_id)

  test_data = data %>%
    filter(year %in% test_years)

  # Allocate an empty vector for cross-validation error
  cv_error = numeric(n.trees)
  for (i in 1:n_folds) {
    # train on the k-1 folds that don't match i, validate on the 1 fold that
    # does (gbm uses the first `train.fraction` proportion as training and the
    # remainder for validation)
    fold_data = rbind(
      train_data[i != train_folds, ],
      train_data[i == train_folds, ]
    )

    if (var(train_data[i != train_folds, "presence"]) > 0) {
      # Fit
      fold_model = gbm(
        presence ~ .,
        distribution = "bernoulli",
        data = fold_data,
        cv.folds = 0,
        train.fraction =  mean(i != train_folds),
        interaction.depth = interaction.depth,
        n.trees = n.trees
      )
      cv_error = cv_error + fold_model$valid.error
    }
  }

  # Find the best number of trees, averaged across the folds
  n_trees_final = which.min(cv_error)

  # Use n.trees to fit a new model to the full training set
  full_model = gbm(
    presence ~ .,
    distribution = "bernoulli",
    data = train_data,
    cv.folds = 0,
    interaction.depth = interaction.depth,
    n.trees = n_trees_final
  )

  p = predict(full_model, newdata = test_data, n.trees = n_trees_final,
              type = "response")

  # Save probabilities and some metadata
  list(species_id = species_id, p = p, n.trees = n_trees_final,
       site_id = test_data$site_id, year = test_data$year)
}


# Run ---------------------------------------------------------------------

sdm_fits = mclapply(unique(occurrences$species_id),
                    fit_species,
                    mc.preschedule = FALSE,
                    mc.cores = mc.cores)
saveRDS(sdm_fits, file = "data/sdm_fits.Rdata")
