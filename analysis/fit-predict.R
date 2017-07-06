# TODO: ensure that sd is correct for all models
# TODO: better names for columns/functions
# TODO: Check namespaces

library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")
dir.create("results", FALSE)

if (!file.exists("observer_model.rds")) {
  fit_observer_model()
}
obs_model = readRDS("observer_model.rds")
observer_variance = mean(obs_model$observer_sigma^2)


# Discard unnecessary data to save memory
bioclim_to_discard = colnames(obs_model$data) %>% 
  grep("^bio", ., value = TRUE) %>% 
  discard(~.x %in% settings$vars) 

# `c()` is needed for a few lines because of distinction between
# vectors and 1D arrays
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect,
         is_future = FALSE) %>% 
  select(-ndvi_ann, -lat, -long, -site_index, -one_of(bioclim_to_discard))

# Reclaim memory
rm(obs_model)
gc()

future = get_env_data(timeframe = "future") %>% 
  filter(site_id %in% !!x_richness$site_id,
         year > !!settings$end_yr) %>% 
  select(which(colnames(.) %in% !!colnames(x_richness))) %>% 
  mutate(richness = 0, iteration = 0,
         observer_effect = 0, is_future = TRUE)

x_richness = bind_rows(x_richness, future)

# Fit & save models --------------------------------------------------------

# "Average" model with observer effects
x_richness %>% 
  filter(!in_train) %>% 
  mutate(mean = intercept + observer_effect + site_effect,
         model = "average", use_obs_model = TRUE) %>% 
  select(site_id, year, mean, sd, iteration, richness, model, use_obs_model) %>% 
  saveRDS(file = "results/avg_TRUE.rds")

# Average model for 2050. Use the observer model to figure out what the
# median observer would have seen at each site (observer effect == 0), then 
# expand to use that same value across all future years.
x_richness %>% 
  mutate(model = "average", use_obs_model = TRUE, 
         mean = intercept + site_effect,
         sd = sqrt(observer_variance + sd^2),
         year = 0, richness = 0) %>% 
  combine_predictions() %>% 
  select(-year, -richness) %>% 
  right_join(select(future, site_id, year), "site_id") %>% 
  saveRDS(file = "results/avg_2050.rds")

# "Average" model without observer effects
# Use site-level means and sds from the training set as test-set predictions
x_richness %>% 
  filter(in_train) %>% 
  group_by(site_id) %>% 
  summarize(mean = mean(richness), sd = sd(richness), model = "average", 
            use_obs_model = FALSE) %>% 
  left_join(select(x_richness, -sd), "site_id") %>% 
  filter(!in_train) %>% 
  select(site_id, year, mean, sd, richness, model, use_obs_model) %>% 
  saveRDS(file = "results/avg_FALSE.rds")


# Forecast-based predictions
# For all combinations of forecast function & use_obs_model (TRUE/FALSE)
# run make_forecasts with data & settings.
expand.grid(fun_name = c("naive", "auto.arima"), 
            use_obs_model = c(TRUE, FALSE),
            stringsAsFactors = FALSE) %>% 
  purrr::transpose() %>% 
  parallel::mclapply(
    function(grid_row){
      do.call(make_all_forecasts, 
              c(x = list(x_richness), settings = list(settings), grid_row))
    },
    mc.cores = 8, 
    mc.preschedule = FALSE
  ) %>% 
  bind_rows() %>% 
  saveRDS(file = "forecast.rds")

# GBM richness with observer effects:
iter_rep(x_richness, make_gbm_predictions, use_obs_model = TRUE) %>% 
  bind_rows() %>% 
  saveRDS(file = "gbm_TRUE.rds")

# GBM richness without observer effects:
iter_rep(x_richness, make_gbm_predictions, use_obs_model = FALSE) %>% 
  bind_rows() %>% 
  saveRDS(file = "gbm_TRUE.rds")

# fit random forest SDMs -------------------------------------------------------

# Get data on individual species occurrences
bbs = get_pop_ts_env_data(settings$start_yr, 
                          settings$end_yr, 
                          settings$last_train_year,
                          settings$min_year_percentage) %>% 
  filter(!is.na(abundance))

# Discard species that don't occur in the training set
bbs = bbs %>% 
  filter(year <= settings$last_train_year) %>% 
  distinct(species_id) %>% 
  left_join(bbs)

dir.create("rf_predictions", showWarnings = FALSE)

gc() # Minimize memory detritus before forking

rf_sdm_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                 settings = settings, use_obs_model = TRUE,
                                 mc.cores = 8) %>% 
  saveRDS(file = "rf_predictions/all_TRUE.rds")

rf_sdm_no_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                    settings = settings, use_obs_model = FALSE,
                                    mc.cores = 8) %>% 
  saveRDS(file = "rf_predictions/all_FALSE.rds")


