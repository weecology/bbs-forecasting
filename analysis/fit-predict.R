# TODO: ensure that sd is correct for all models

library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")

obs_model = readRDS("observer_model.rds")

iters = 1:25 # During interactive use/testing, only include a few MCMC samples

# `c()` is needed for a few lines because of distinction between 
# vectors and 1D arrays, I think?
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect) %>% 
  filter(iteration %in% iters)


average_model_predictions = x_richness %>% 
  filter(!in_train) %>% 
  mutate(mean = intercept + observer_effect + site_effect,
         model = "average", obs_model = TRUE) %>% 
  select(site_id, year, mean, sd, iteration, richness, model, obs_model)

average_no_obs = x_richness %>% 
  filter(in_train) %>% 
  group_by(site_id) %>% 
  summarize(mean = mean(richness), sd = sd(richness), model = "average", obs_model = FALSE) %>% 
  left_join(select(x_richness, -sd), "site_id") %>% 
  filter(!in_train) %>% 
  select(site_id, year, mean, sd, richness, model, obs_model)

naive_model_predictions = make_all_forecasts(x_richness, 
                                             "naive", 
                                             obs_model = TRUE, 
                                             settings = settings)
naive_no_obs = make_all_forecasts(x_richness, 
                                  "naive", 
                                  obs_model = FALSE, 
                                  settings = settings)

auto_model_predictions = make_all_forecasts(x_richness, 
                                            "auto.arima", 
                                            obs_model = TRUE, 
                                            settings = settings,
                                            seasonal = FALSE)

auto_no_obs = make_all_forecasts(x_richness, "auto.arima", 
                                 obs_model = FALSE, 
                                 settings = settings,
                                 seasonal = FALSE)

gbm_richness_predictions = x_richness %>%
  group_by(iteration) %>% 
  by_slice(make_gbm_predictions, obs_model = TRUE, .collate = "rows")

# Don't need to group/by_slice because only fitting one iteration
gbm_no_obs = x_richness %>%
  filter(iteration == 1) %>%
  make_gbm_predictions(obs_model = FALSE)

p = bind_rows(average_model_predictions, average_no_obs, 
              naive_model_predictions, naive_no_obs,
              auto_no_obs, auto_model_predictions,
              gbm_richness_predictions, gbm_no_obs)

p %>% 
  filter(!is.na(richness)) %>% 
  combine_predictions() %>% 
  mutate(`mean deviance` = -2 * dnorm(richness, mean, sd, log = TRUE),
         mae = abs(richness - mean), 
         obs_model = forcats::fct_relevel(factor(obs_model), "TRUE")) %>% 
  ggplot(aes(x = year, y = `mean deviance`, color = model, linetype = obs_model)) + 
  geom_smooth(method = "gam", formula = y ~ s(x)) + 
  cowplot::theme_cowplot() + 
  scale_color_brewer(type = "qual", palette = "Dark2")
