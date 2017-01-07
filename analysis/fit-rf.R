library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")

if (!file.exists("observer_model.rds")) {
  fit_observer_model()
} else {
  obs_model = readRDS("observer_model.rds")
}

bbs = get_bbs_data()

# `c()` is needed for a few lines because of distinction between 
# vectors and 1D arrays, I think?
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect)

rm(obs_model)
gc()


rf_sdm_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                 settings = settings, obs_model = TRUE)

saveRDS(rf_sdm_obs, file = "rf_predictions/all_TRUE.rds")
rm(rf_sdm_obs)
gc()

rf_sdm_no_obs = rf_predict_richness(bbs = bbs, x_richness = x_richness, 
                                    settings = settings, obs_model = FALSE) 

saveRDS(rf_sdm_no_obs, file = "rf_predictions/all_FALSE.rds")

