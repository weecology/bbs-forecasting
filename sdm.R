# From forecasting-bbs-R.ipynb --------------------------------------------
library(forecast)
library(Hmisc)
library(dplyr)
library(broom)
library(ggplot2)
library(tidyr)
library(mgcv)
source("forecast-bbs-core.R")
library(sp)
library(raster)
library(maptools)
library(rgeos)
library(rgdal)

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25


# start SDM ---------------------------------------------------------------

n_test_years = 5
interaction.depth = 8
n.trees = 1E4

train_years = start_yr:(end_yr - n_test_years)
test_years = (end_yr - n_test_years + 1):end_yr

library(gbm)
library(viridis)

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
  occurrences %>%
    filter(species_id == x) %>%
    dplyr::select(site_id, year, abundance) %>%
    right_join(sampling_events, by = c("site_id", "year")) %>%
    replace_na(list(abundance = 0))
}



# Generate fold_id --------------------------------------------------------

set.seed(0)

locations = sampling_events %>%
  distinct(site_id) %>%
  dplyr::select(site_id, lat, long)

cluster_id = kmeans(locations[ , 2:3], n_clusters, nstart = 100)$cluster

fold_matrix = matrix(sample(1:n_clusters), ncol = n_folds)

fold_id = numeric(nrow(locations))
for (i in 1:n_folds) {
  fold_id = fold_id + i * (cluster_id %in% fold_matrix[, i])
}

plot(
  locations[,3:2],
  col = fold_id,
  pch = 16
)

locations = cbind(locations, fold_id = fold_id)

# fit gbm -----------------------------------------------------------------
fit_species = function(species_id){

  data = make_species_data(x = species_id) %>%
    mutate(presence = abundance > 0)


  stopifnot(all(data$lat == sampling_events$lat))

  train_data = data %>%
    filter(year %in% train_years) %>%
    inner_join(locations[ , c("fold_id", "site_id")], by = "site_id") %>%
    dplyr::select(-lat, -long, -site_id, -year, -abundance)

  train_folds = train_data$fold_id

  train_data = dplyr::select(train_data, -fold_id)

  test_data = data %>%
    filter(year %in% test_years)

  cv_error = numeric(n_folds)
  for (i in 1:n_folds) {
    print(i)
    # train on the k-1 folds that don't match i, validate on the 1 fold that does
    fold_data = rbind(
      train_data[i != train_folds, ],
      train_data[i == train_folds, ]
    )
    fold_model = gbm(
      presence ~ .,
      distribution = "bernoulli",
      data = fold_data,
      cv.folds = 0,
      train.fraction =  mean(i != train_folds),
      interaction.depth = interaction.depth,
      n.trees = n.trees
    )
    cv_error = cv_error + fold_model$valid.error / n_folds
  }

  n.trees = which.min(cv_error)

  full_model = gbm(
    presence ~ .,
    distribution = "bernoulli",
    data = train_data,
    cv.folds = 0,
    interaction.depth = interaction.depth,
    n.trees = n.trees
  )

  p = predict(full_model, newdata = test_data, n.trees = n.trees,
              type = "response")

  list(p = p, n.trees = n.trees)
}

