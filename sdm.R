# From forecasting-bbs-R.ipynb --------------------------------------------
devtools::load_all()
start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25

# start SDM ---------------------------------------------------------------
library(gbm)
library(parallel)
library(dplyr)
library(broom)
library(tidyr)

set.seed(0)
mc.cores = 2
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
occurrences = get_bbs_data(start_yr = 1982, end_yr = 2013, min_num_yrs = 25)

env = get_env_data() %>%
  filter_ts(start_yr, end_yr, min_num_yrs) %>%
  inner_join(distinct(occurrences, site_id, lat, long, year),
             by = c("site_id", "year"))

# Determine presence-absence ----------------------------------------------

x = tidyr::expand(occurrences, site_id, year) %>%
  inner_join(env, by = c("site_id", "year")) %>%
  arrange(year, site_id)

locations = distinct(x, site_id, lat, long)

y = distinct(occurrences, year, species_id, site_id, .keep_all = TRUE) %>%
  inner_join(select(x, year, site_id), c("year", "site_id")) %>%
  mutate(presence = abundance > 0) %>%
  select(-abundance, lat, long, start_time, month, day) %>%
  spread(key = species_id, value = presence, fill = 0) %>%
  arrange(year, site_id)


# Spatial cross validation ------------------------------------------------
# assign many clusters of routes that are close in lat/long
cluster_id = kmeans(locations[ , 2:3], n_clusters, nstart = 100)$cluster

# Randomly combine the clusters to make the cross-validation folds
fold_matrix = matrix(sample(1:n_clusters), ncol = n_folds)
fold_id = numeric(nrow(locations))
for (i in 1:n_folds) {
  fold_id = fold_id + i * (cluster_id %in% fold_matrix[, i])
}

x = x %>%
  inner_join(
    data_frame(site_id = locations$site_id, fold_id = fold_id),
    "site_id"
  )

stopifnot(x$site_id == y$site_id, x$year == y$year)

train_x = filter(x, year %in% train_years)
test_x = filter(x, !(year %in% train_years))
train_y = filter(y, year %in% train_years)
test_y = filter(y, !(year %in% train_years))

fold_ids = train_x$fold_id  #plot(locations[,3:2], col = fold_ids, pch = fold_id)

cv_train_rows = function(i){
  i != fold_ids
}
cv_test_rows = function(i){
  i == fold_ids
}
# fit gbm -----------------------------------------------------------------

# Function to cross-validate a GBM model with different numbers of trees, then
# fit the optimal model size to the full training set.
fit_species = function(species_id){
  species_id = as.character(species_id)
  train_data = cbind(train_x, y = train_y[[species_id]]) %>%
    select(-site_id, -year, -lat, -long, -fold_id)

  # Allocate an empty vector for cross-validation error
  cv_error = numeric(n.trees)
  for (i in 1:n_folds) {
    fold_data = rbind(
      train_data[cv_train_rows(i), ],
      train_data[cv_test_rows(i), ]
    )

    if (var(train_data[cv_train_rows(i), "y"]) > 0) {
      fold_model = gbm(
        y ~ .,
        distribution = "bernoulli",
        data = fold_data,
        cv.folds = 0,
        train.fraction =  mean(cv_train_rows(i)),
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
    y ~ .,
    distribution = "bernoulli",
    data = train_data,
    cv.folds = 0,
    interaction.depth = interaction.depth,
    n.trees = n_trees_final
  )

  p = predict(full_model, newdata = test_x, n.trees = n_trees_final,
              type = "response")

  # Save probabilities and some metadata
  data.frame(species_id = species_id, p = p, n.trees = n_trees_final,
             site_id = test_x$site_id, year = test_x$year)
}


# Run ---------------------------------------------------------------------

sdm_fits = mclapply(
  unique(occurrences$species_id)[1:2],
  fit_species,
  mc.preschedule = FALSE,
  mc.cores = mc.cores
) %>%
  bind_rows()
saveRDS(sdm_fits, file = "data/sdm_fits.Rdata")
