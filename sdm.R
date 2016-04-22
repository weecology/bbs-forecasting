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



# Generate fold_id --------------------------------------------------------

set.seed(0)

locations = env %>% dplyr::select(site_id, lat, long) %>% distinct()

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

# fit gbm -----------------------------------------------------------------

# NOTE: NOT USING THE CV FOLDS DEFINED ABOVE!

data = make_species_data(x = 2890) %>%
  mutate(presence = abundance > 0)
g = gbm(
  presence ~ .,
  distribution = "bernoulli",
  data = data %>% dplyr::select(-lat, -long, -site_id, -year, -abundance),
  cv.folds = n_folds,
  interaction.depth = 4,
  n.trees = 1E3
)

p = predict(g, type = "response", n.trees = g$n.trees)
plot(
  jitter(data$lat, amount = 1/2) ~ jitter(data$long, amount = 1/2), 
  col = rep(viridis(length(p) / 100), each = 100)[order(p)],
  pch = ".",
  cex = 8
)

results = data %>% 
  cbind(p = p) %>% 
  dplyr::group_by(site_id, lat, long) %>% 
  summarise(count = sum(presence), p = mean(p)) 

ggplot(results, aes(x = long, y = lat, color = count)) + geom_point(size = 4) + scale_color_viridis()
ggplot(results, aes(x = long, y = lat, color = p)) + geom_point(size = 4) + scale_color_viridis()
