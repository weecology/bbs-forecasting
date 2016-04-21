library(magrittr)
# Parameters for spatial cross-validation
n_clusters = 100
n_folds = 10
stopifnot(n_clusters %% n_folds == 0)


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
  left_join(env, by = c("site_id", "year"))

x = 2890

abundance = occurrences %>% 
  filter(species_id == x) %>% 
  dplyr::select(site_id, year, abundance) %>%
  right_join(sampling_events, by = c("site_id", "year")) %>%
  replace_na(list(abundance = 0))

