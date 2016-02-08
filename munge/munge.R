library(raster)       # for pointDistance
library(dplyr)        # for data manipulation
library(magrittr)     # for data manipulation
library(tidyr)        # for data manipulation
library(geosphere)    # for regularCoordinates


focal_year = 1990     # use data through this year to predict the rest
min_occur = 25        # include a species if it was observed at least `min_occur` 
                      # times by focal_year


database = src_sqlite("data/bbsforecasting.sqlite")

bbs = tbl(database, "bbs_data") %>% collect()


# train-test split --------------------------------------------------------

# All the points within inner_radius of a center point will be in the test set.
# Everything more than outer_radius away will be in the training set.
# radii are in meters
centers = regularCoordinates(16)
inner_radius = 1.0E5
outer_radius = 2.5E5

site_lat_long = bbs %>% dplyr::select(site_id, lat, long) %>% distinct()

dists = pointDistance(
  centers, 
  cbind(site_lat_long$long, site_lat_long$lat), 
  longlat = TRUE
)

site_lat_long$in_train = apply(dists, 2, min) > outer_radius
site_lat_long$in_test = apply(dists, 2, min) < inner_radius

bbs = inner_join(bbs, site_lat_long)


# Drop the rarest species -------------------------------------------------

# Only include a species if it was obseved along `min_occur` or more runs
# by `focal_year`.  Otherwise, there won't be enough data to make predictions
included_species = bbs %>% 
  filter(year <= focal_year) %>% 
  group_by(species_id) %>% 
  summarise(include = n() >= min_occur) %>% 
  filter(include) %>%
  extract2("species_id")


bbs = bbs %>% filter(species_id %in% included_species)


# prism -------------------------------------------------------------------

site_year = bbs[ , c("site_id", "year")] %>% distinct

prism =
  tbl(database, 'prism_bbs_data') %>%
  inner_join(site_year, by = c("site_id", "year"), copy = TRUE) %>%
  collect() %>% 
  mutate(clim_var_month = paste(clim_var, month, sep="_")) %>%
  select(-clim_var, -month) %>%
  spread(clim_var_month, value)

# output ------------------------------------------------------------------

full_bbs = bbs %>% 
  spread(species_id, abundance, fill = 0) %>% 
  inner_join(prism, by = c("site_id", "lat", "long", "year"))



saveRDS(full_bbs, file = "data/munged.rds")
saveRDS(focal_year, "data/focal_year.rds")
