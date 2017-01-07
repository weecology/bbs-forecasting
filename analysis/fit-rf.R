library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")

if (!file.exists("observer_model.rds")) {
  fit_observer_model()
} else {
  obs_model = readRDS("observer_model.rds")
}

rf_sdm_obs = rf_predict_richness(bbs = get_bbs_data(), x_richness = x_richness, 
                                 settings = settings, obs_model = TRUE)

saveRDS(rf_sdm_obs, file = "rf_predictions/all_TRUE.rds")
rm(rf_sdm_obs)
gc()

rf_sdm_no_obs = rf_predict_richness(bbs = get_bbs_data(), x_richness = x_richness, 
                                    settings = settings, obs_model = FALSE) 

saveRDS(rf_sdm_no_obs, file = "rf_predictions/all_FALSE.rds")

