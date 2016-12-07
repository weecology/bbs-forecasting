combine_predictions = function(x){
  x %>% 
    group_by(site_id, year, model, obs_model, richness) %>% 
    summarize(sd = sqrt(mean(var(mean) + mean(sd^2))), 
              mean = mean(mean)) %>% 
    ungroup()
}


# TODO: ensure that sd is correct for all models

library(tidyverse)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")

obs_model = readRDS("observer_model.rds")

iters = 1:50 # During interactive use/testing, only include a few MCMC samples

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


make_naive_predictions = function(x){
  # Set the `level` so that `upper` and `lower` are 2 sd apart.
  fc = forecast::naive(x$expected_richness, h = settings$end_yr - settings$last_train_year,
                  level = pnorm(0.5))
  # Distance between `upper` and `lower` is 2 sd, so divide by 2
  tibble(year = seq(settings$last_train_year + 1, settings$end_yr), 
         mean = c(fc$mean), sd = c(fc$upper - fc$lower) / 2, model = "naive",
         obs_model = TRUE)
}

naive_model_predictions = x_richness %>% 
  filter(year <= settings$last_train_year) %>% 
  complete(site_id, year) %>% 
  group_by(site_id, iteration) %>% 
  arrange(year) %>% 
  by_slice(make_naive_predictions, .collate = "row") %>% 
  left_join(select(x_richness, -sd), c("site_id", "year", "iteration")) %>% 
  mutate(mean = mean + observer_effect) %>% 
  select(site_id, year, mean, sd, iteration, richness, model, obs_model)


make_naive_predictions2 = function(x){
  # Set the `level` so that `upper` and `lower` are 2 sd apart.
  fc = forecast::naive(x$richness, h = settings$end_yr - settings$last_train_year,
                       level = pnorm(0.5))
  # Distance between `upper` and `lower` is 2 sd, so divide by 2
  tibble(year = seq(settings$last_train_year + 1, settings$end_yr), 
         mean = c(fc$mean), sd = c(fc$upper - fc$lower) / 2, model = "naive",
         obs_model = FALSE)
}
naive_no_obs = x_richness %>% 
  filter(year <= settings$last_train_year) %>% 
  complete(site_id, year) %>% 
  group_by(site_id, iteration) %>% 
  arrange(year) %>% 
  by_slice(make_naive_predictions2, .collate = "row") %>% 
  left_join(select(x_richness, -sd), c("site_id", "year", "iteration")) %>% 
  mutate(mean = mean) %>% 
  select(site_id, year, mean, sd, iteration, richness, model, obs_model)




make_gbm_predictions = function(x){
  train = filter(x, year <= settings$last_train_year)
  g = gbm::gbm(expected_richness ~ 
                 bio2 + bio3 + bio5 + bio8 + bio9 + bio15 + bio16 + bio18 + 
                 ndvi_sum + 
                 ndvi_win + 
                 elevs, 
               data = train,
               distribution = "gaussian",
               interaction.depth = 3,
               shrinkage = .01,
               n.trees = 5000)
  test = filter(x, year > settings$last_train_year)
  mean = predict(g, test, n.trees = gbm::gbm.perf(g, plot.it = FALSE)) + 
    test$observer_effect
  cbind(test, mean = mean, model = "richness_gbm",  obs_model = TRUE, 
        stringsAsFactors = FALSE) %>% 
    select(site_id, year, mean, sd, richness, model, obs_model)
}

# gbm_predictions = x_richness %>% 
#   group_by(iteration) %>% 
#   by_slice(make_gbm_predictions, .collate = "row")

p = bind_rows(average_model_predictions, average_no_obs, naive_model_predictions, naive_no_obs)

square = function(x)x^2

p %>% 
  filter(!is.na(richness)) %>% 
  combine_predictions() %>% 
  mutate(`mean deviance` = -2 * dnorm(richness, mean, sd, log = TRUE),
         mae = abs(richness - mean), 
         obs_model = forcats::fct_relevel(factor(obs_model), "TRUE")) %>% 
  ggplot(aes(x = year, y = `mean deviance`, color = model, linetype = obs_model)) + 
  geom_smooth(method = "gam") + 
  cowplot::theme_cowplot()
