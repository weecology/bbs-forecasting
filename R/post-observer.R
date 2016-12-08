make_forecast = function(x, fun_name, obs_model, settings, ...){
  # `Forecast`` functions want NAs for missing years, & want the years in order
  x = complete(x, year = settings$start_yr:settings$last_train_year) %>% 
    arrange(year)
  
  if (obs_model) {
    # observer effect will be added back in `make_all_forecasts`
    response_variable = "expected_richness"
  } else {
    response_variable = "richness"
  }
  
  h = settings$end_yr - settings$last_train_year
  
  # Set the `level` so that `upper` and `lower` are 2 sd apart.
  level = pnorm(0.5)
  
  if (all(is.na(x[[response_variable]]))) {
    warning("Empty response variable")
    return(list(NA))
  }
  
  
  if (fun_name == "naive") {
    # `forecast::naive` refuses to predict when the final observation is NA.
    # But we can fit the same model with `Arima(order = c(0,1,0))`
    fun = partial(Arima, order = c(0, 1, 0))
  } else {
    # Just get the named function
    fun = getFromNamespace(fun_name, "forecast")
  }
  
  fc = fun(y = x[[response_variable]], ...) %>% 
    forecast::forecast(h = h, level = level)
  
  if (any(is.na(fc$mean))) {
    browser()
  }
  
  # Distance between `upper` and `lower` is 2 sd, so divide by 2
  tibble(year = seq(settings$last_train_year + 1, settings$end_yr), 
         mean = c(fc$mean), sd = c(fc$upper - fc$lower) / 2, model = fun_name,
         obs_model = obs_model)
}

make_all_forecasts = function(x, fun_name, obs_model, 
                              settings, ...){
  forecast_data = x %>% 
    filter(year <= settings$last_train_year) %>% 
    group_by(site_id, iteration)
  
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

make_gbm_predictions = function(x, obs_model){
  train = filter(x, year <= settings$last_train_year)
  test = filter(x, year > settings$last_train_year)
  
  if (obs_model) {
    train$y = train$expected_richness
  } else {
    train$y = train$richness
  }
  
  g = gbm::gbm(y ~ 
                 bio2 + bio3 + bio5 + bio8 + bio9 + bio15 + bio16 + bio18 + 
                 ndvi_sum + 
                 ndvi_win + 
                 elevs, 
               data = train,
               distribution = "gaussian",
               interaction.depth = 3,
               shrinkage = .01,
               n.trees = 5000)
  
  mean = predict(g, test, n.trees = gbm::gbm.perf(g, plot.it = FALSE))
  if (obs_model) {
    mean = mean + test$observer_effect
  }
  
  cbind(test, mean = mean, model = "richness_gbm",  obs_model = obs_model, 
        stringsAsFactors = FALSE) %>% 
    select(site_id, year, mean, sd, richness, model, obs_model)
}


combine_predictions = function(x){
  # If we only have one iteration, we want to say the variance is zero, not that
  # it's NA.
  safe_var = function(x){
    if (length(x) == 1) {
      0
    } else {
      var(x)
    }
  }
  
  x %>% 
    group_by(site_id, year, model, obs_model, richness) %>% 
    summarize(sd = sqrt(mean(safe_var(mean) + mean(sd^2))), 
              mean = mean(mean)) %>% 
    ungroup()
}
