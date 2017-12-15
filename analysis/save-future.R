devtools::load_all()
library(tidyverse)
library(cowplot)
library(ggjoy)
library(mvtnorm)

settings = yaml::yaml.load_file("settings.yaml")

timeframe = "train_32"
prepend_timeframe = function(x) {
  paste0("results", "/", timeframe, "/", x)
}


forecast_results = readRDS(prepend_timeframe("forecast.rds"))
gbm_results = bind_rows(readRDS(prepend_timeframe("gbm_TRUE.rds")))
mistnet_results = lapply(dir(prepend_timeframe("mistnet_output/"),
                             full.names = TRUE,
                             pattern = "TRUE"),
                         readRDS) %>%
  bind_rows()


average_results = bind_rows(readRDS(prepend_timeframe("avg_TRUE.rds"))) %>% 
  mutate(model = "average", use_obs_model = TRUE)

my_var = function(x){
  ifelse(length(x) == 1, 0, var(x))
}


rf_results = bind_rows(readRDS(prepend_timeframe("rf_predictions/all_TRUE.rds")))

bound = bind_rows(
  forecast_results, 
  gbm_results, 
  average_results, 
  rf_results, 
  mistnet_results
) %>% 
  group_by(site_id, year, model, use_obs_model) %>% 
  summarize(sd = sqrt(mean(sd^2 + my_var(mean))), mean = mean(mean)) %>% 
  ungroup() %>% 
  mutate(
    model_name = forcats::fct_recode(
      factor(model),
      Average = "average",
      Naive = "naive",
      Auto.arima = "auto.arima",
      `Stacked random forest SDMs` = "rf_sdm",
      `Mistnet JSDM` = "mistnet",
      `GBM richness regression` = "richness_gbm"
    )
  )

bound %>% 
  filter(use_obs_model) %>% 
  filter(year > settings$timeframes$train_32$last_train_year) %>% 
  select(-use_obs_model) %>% 
  write_csv("results/train_32/summaries.csv.gz")
