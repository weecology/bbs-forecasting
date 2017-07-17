library(tidyverse)
library(purrrlyr) # for by_slice
library(parallel)
devtools::load_all()
timeframe = "train_32"


# Ingest the `settings` and put the timeframe of interest at the top level
# of the list.
settings = yaml::yaml.load_file("settings.yaml")
settings = c(settings, settings$timeframes[[timeframe]])
settings$timeframes = NULL


prepend_timeframe = function(x) {
  paste0("results", "/", timeframe, "/", x)
}
dir.create(prepend_timeframe(NULL), showWarnings = FALSE, recursive = TRUE)

if (!file.exists(prepend_timeframe("observer_model.rds"))) {
  fit_observer_model(settings = settings, output_file = prepend_timeframe("observer_model.rds"))
}
obs_model = readRDS(prepend_timeframe("observer_model.rds"))

# Discard unnecessary data to save memory
bioclim_to_discard = colnames(obs_model$data) %>% 
  grep("^bio", ., value = TRUE) %>% 
  discard(~.x %in% settings$vars) 

# `c()` is needed for a few lines because of distinction between
# vectors and 1D arrays
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect) %>% 
  select(-ndvi_ann, -lat, -long, -site_index, -one_of(bioclim_to_discard))


if (settings$timeframe == "future") {
  future = get_env_data(timeframe = "future") %>% 
    filter(site_id %in% !!x_richness$site_id,
           year > !!settings$last_train_year) %>% 
    select(which(colnames(.) %in% !!colnames(x_richness)))
  
  # Standard deviation of "future" observers for each iteration, add as a 
  # column in x_richness
  observer_sigmas = obs_model$observer_sigma
  x_richness = enframe(observer_sigmas, name = "iteration", 
                       value = "observer_sigma") %>% 
    right_join(x_richness, "iteration")
} else {
    future = NULL
    observer_sigmas = NULL
}

# Reclaim memory
rm(obs_model)
gc()


# Fit & save models --------------------------------------------------------



if (settings$timeframe == "future") {
  year_grid = expand.grid(iteration = unique(x_richness$iteration), 
                          year = years_to_use)
  
  # Calculate expected richness for "typical" observer, add in noise for
  # unknown observers. Then repeat for all year/iteration combinations
  x_richness %>% 
    mutate(mean = intercept + site_effect,
           sd = sqrt(sd^2 + observer_sigma^2)) %>% 
    distinct(iteration, site_id, mean, sd) %>% 
    left_join(year_grid, "iteration") %>% 
    saveRDS(file = prepend_timeframe("avg_TRUE.rds"))
  
  # When there's no observer effect, just calculate the empirical mean and sd 
  # for each site, then replicate across years.
  avg_false = x_richness %>% 
    filter(in_train) %>% 
    group_by(site_id) %>% 
    summarize(mean = mean(richness), sd = sd(richness))
  map_df(years_to_use, function(x) cbind(avg_false, year = x)) %>% 
    saveRDS(file = prepend_timeframe("avg_FALSE.rds"))
} else {
  # "Average" model with observer effects
  x_richness %>% 
    filter(!in_train) %>% 
    mutate(mean = intercept + observer_effect + site_effect,
           model = "average", use_obs_model = TRUE) %>% 
    select(site_id, year, mean, sd, iteration, richness, model, use_obs_model) %>% 
    saveRDS(file = prepend_timeframe("avg_TRUE.rds"))
  
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
    saveRDS(file = prepend_timeframe("avg_FALSE.rds"))
}


# Forecast-based predictions
# For all combinations of forecast function & use_obs_model (TRUE/FALSE)
# run make_forecasts with data & settings.
expand.grid(fun_name = c("naive", "auto.arima"), 
            use_obs_model = c(TRUE, FALSE),
            stringsAsFactors = FALSE) %>% 
  transpose() %>% 
  mclapply(
    function(grid_row){
      make_all_forecasts(x = x_richness, settings = settings, 
                         observer_sigmas = observer_sigmas,
                         fun_name = grid_row$fun_name, 
                         use_obs_model = grid_row$use_obs_model)
    },
    mc.cores = 8, 
    mc.preschedule = FALSE
  ) %>% 
  bind_rows() %>% 
  saveRDS(file = prepend_timeframe("forecast.rds"))

# GBM richness with observer effects:
x_richness %>%
  split(., .$iteration) %>% 
  mclapply(make_gbm_predictions, 
           use_obs_model = TRUE,
           settings = settings, 
           future = future,
           observer_sigmas = NULL,
           mc.cores = 8,
           mc.preschedule = FALSE) %>% 
  bind_rows() %>% 
  saveRDS(file = prepend_timeframe("gbm_TRUE.rds"))

# GBM richness without observer effects:
x_richness %>%
  split(., .$iteration) %>% 
  mclapply(make_gbm_predictions, 
           use_obs_model = FALSE,
           settings = settings, 
           future = future,
           observer_sigmas = NULL,
           mc.cores = 8,
           mc.preschedule = FALSE) %>% 
  bind_rows() %>% 
  saveRDS(file = prepend_timeframe("gbm_FALSE.rds"))

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

rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                    settings = settings, use_obs_model = TRUE,
                    future = future, 
                    observer_sigmas = observer_sigmas) %>% 
  saveRDS(file = prepend_timeframe("rf_predictions/all_TRUE.rds"))

rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                    settings = settings, use_obs_model = FALSE,
                    future = future, 
                    observer_sigmas = observer_sigmas) %>% 
  saveRDS(file = prepend_timeframe("rf_predictions/all_FALSE.rds"))
