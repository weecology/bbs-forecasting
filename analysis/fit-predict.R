# If we only have one iteration, we want to say the variance is zero, not that
# it's NA.
safe_var = function(x){
  if (length(x) == 1) {
    0
  } else {
    var(x)
  }
}

combine_predictions = function(x){
  x %>% 
    group_by(site_id, year, model, obs_model, richness) %>% 
    summarize(sd = sqrt(mean(safe_var(mean) + mean(sd^2))), 
              mean = mean(mean)) %>% 
    ungroup()
}


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


make_forecast = function(x, fun_name, obs_model, settings, ...){
  if (obs_model) {
    # observer effect will be added back in `make_all_forecasts`
    response_variable = "expected_richness"
  } else {
    response_variable = "richness"
  }
  fun = getFromNamespace(fun_name, "forecast")
  
  # Set the `level` so that `upper` and `lower` are 2 sd apart.
  fc = fun(x[[response_variable]], ...) %>% 
    forecast::forecast(h = settings$end_yr - settings$last_train_year,
                       level = pnorm(0.5))
  # Distance between `upper` and `lower` is 2 sd, so divide by 2
  tibble(year = seq(settings$last_train_year + 1, settings$end_yr), 
         mean = c(fc$mean), sd = c(fc$upper - fc$lower) / 2, model = fun_name,
         obs_model = obs_model)
}

make_all_forecasts = function(x, fun_name, obs_model, 
                              settings, ...){
  # `Forecast`` functions want NAs for missing years, and want the years in order
  forecast_data = x %>% 
    filter(year <= settings$last_train_year) %>% 
    complete(site_id, year) %>%
    group_by(site_id, iteration) %>%
    arrange(year)
  
  if (!obs_model) {
    # Without an observation model, all the iterations will be the same.
    # Don't bother fitting the same model to each iteration
    forecast_data = filter(forecast_data, iteration == 1)
  }
  
  # TODO: by_slice will be deprecated. See if I can use mutate_all instead?
  # Dropping x_richness$sd before joining so it can be replaced by forecast sd.
  out = by_slice(forecast_data, make_forecast, fun_name = fun_name, 
           obs_model = obs_model, settings = settings, ...,
           .collate = "row") %>%
    left_join(select(x_richness, -sd), c("site_id", "year", "iteration"))
  
  if (obs_model) {
    # Observer effect was subtraced out in make_forcast. Add it back in here.
    out = mutate(out, mean = mean + observer_effect)
  }

  select(out, site_id, year, mean, sd, iteration, richness, model, obs_model)
}

naive_model_predictions = make_all_forecasts(x_richness, "naive", 
                                             obs_model = TRUE, 
                                             settings = settings)
naive_no_obs = make_all_forecasts(x_richness, "naive", 
                                             obs_model = FALSE, 
                                             settings = settings)

auto_no_obs = make_all_forecasts(x_richness, "auto.arima", 
                                 obs_model = FALSE, 
                                 settings = settings,
                                 seasonal = FALSE)


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

p = bind_rows(average_model_predictions, average_no_obs, 
              naive_model_predictions, naive_no_obs,
              auto_no_obs)

p %>% 
  filter(!is.na(richness)) %>% 
  combine_predictions() %>% 
  mutate(`mean deviance` = -2 * dnorm(richness, mean, sd, log = TRUE),
         mae = abs(richness - mean), 
         obs_model = forcats::fct_relevel(factor(obs_model), "TRUE")) %>% 
  ggplot(aes(x = year, y = `mean deviance`, color = model, linetype = obs_model)) + 
  geom_smooth(method = "gam") + 
  cowplot::theme_cowplot()
