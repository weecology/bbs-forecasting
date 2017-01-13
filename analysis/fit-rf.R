set.seed(1)
library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")

if (!file.exists("observer_model.rds")) {
  obs_model = fit_observer_model()
} else {
  obs_model = readRDS("observer_model.rds")
}

bbs = get_pop_ts_env_data(settings$start_yr, 
                          settings$end_yr, 
                          settings$min_num_yrs) %>% 
  filter(!is.na(abundance))

# Discard unnecessary data to save memory
bioclim_to_discard = colnames(obs_model$data) %>% 
  grep("^bio", ., value = TRUE) %>% 
  discard(~.x %in% settings$vars)

# `c()` is needed for a few lines because of distinction between 
# vectors and 1D arrays, I think?
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect) %>% 
  select(-ndvi_ann, -lat, -long, -site_index, -expected_richness, -sd,
         -intercept, -site_effect, -one_of(bioclim_to_discard))

rm(obs_model)
gc()


rf_sdm_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                 settings = settings, use_obs_model = TRUE,
                                 mc.cores = 8)

dir.create("rf_predictions", showWarnings = FALSE)
saveRDS(rf_sdm_obs, file = "rf_predictions/all_TRUE.rds")
rm(rf_sdm_obs)
gc()

rf_sdm_no_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                    settings = settings, use_obs_model = FALSE,
                                    mc.cores = 8) 

saveRDS(rf_sdm_no_obs, file = "rf_predictions/all_FALSE.rds")

