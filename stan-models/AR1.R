library(dplyr)
library(tidyr)
library(rstan)
devtools::load_all()

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25
last_training_year <- 2003


raw = get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
  select(site_id, year, richness, observer_id) %>%
  complete(site_id, year) %>%
  arrange(site_id, year)

training_observations = filter(raw, year <= last_training_year)
testing_observations = filter(raw, year > last_training_year)

training_observations %>%
  group_by(year, site_id) %>%
  tidyr::spread(key = year, value = richness)

data = training_observations %>%
  mutate(observation_index = 1:nrow(.)) %>%
  rename(observed_richness = richness) %>%
  filter(!is.na(observed_richness)) %>%
  as.list()

center = attr(scale(data$observed_richness), "scaled:center")
scale = attr(scale(data$observed_richness), "scaled:scale")

# Observation-level pointers
# Note that site-level random effect includes non-observed years, so its index
# must as well
data$observed_richness = c(scale(data$observed_richness))
data$observer_index = as.integer(factor(data$observer_id))
data$site_index = as.integer(factor(training_observations$site_id))

# Vector lengths
data$N_observations = length(data$observed)
data$N_train_years = length(start_yr:last_training_year)
data$N_test_years = length(seq(last_training_year + 1, end_yr))
data$N_sites = length(unique(data$site_id))
data$N_observers = max(data$observer_index)

# year pointers (based on _all_ site/year combinations, not just ones with
# observations)
data$which_first = which(training_observations$year == start_yr)
data$which_non_first = which(training_observations$year != start_yr)

# Drop unused variables
data$site_id = NULL
data$year = NULL
data$observer_id = NULL

m = stan_model(
  "stan-models/AR1.stan"
)

# Check the number of physical cores on the machine. With 2 or fewer cores, just
# use one core.  Otherwise, use all of them.
cores = ifelse(parallel::detectCores(logical = FALSE) <= 2, 1, parallel::detectCores(logical = FALSE))

model = sampling(m, data = data, control = list(adapt_delta = 0.9),
                 cores = cores, chains = 1, refresh = 1)
saveRDS(model, file = "stan-models/model-object.RDS")

# extract model estimates -------------------------------------------------
unscale = function(data){
  data * scale + center
}
get_quantiles = function(x, q = c(.025, .975)){
  t(apply(x, 2, quantile, q))
}

extracted = rstan::extract(model)

expected = structure(extracted$y, dim = c(nrow(extracted$y), data$N_train_years, data$N_sites))
#future_y = array(extracted$future_y, dim = c(nrow(extracted$future_y), data$N2, length(data$site_length)))

k = sample.int(dim(expected)[3], 1)
matplot(start_yr:last_training_year, unscale(t(expected[,,k])), col = "#00000010", type = "l", lty = 1, ylab = "richness", xlab = "year", xlim = c(start_yr, end_yr), ylim = range(unscale(data$observed)))
#matlines(seq(last_training_year+1, end_yr), unscale(t(future_y[,,k])), col = "#00000010", type = "l", lty = 1, ylab = "richness", xlab = "year")
point_data = training_observations %>% filter(as.integer(factor(site_id)) == k) %>%
  mutate(observer_id = factor(observer_id))
point_data %>% lines(richness ~ year, data = ., col = 2, lwd = 2)
point_data %>%
  points(richness ~ year, data = ., col = 2 + as.integer(observer_id), pch = 16, ylim = range(unscale(data$observed)))


# Vestigial pre-observer code ---------------------------------------------
stop()
# evaluation --------------------------------------------------------------

diffs = as.matrix(predicted_means - testing_observations)

# RMSE
sqrt(mean(diffs^2, na.rm = TRUE))

# Proportion outside the confidence interval (should each be close to 0.025)
mean(testing_observations < lower_CIs, na.rm = TRUE) # below lower
mean(testing_observations > upper_CIs, na.rm = TRUE) # above upper

# Proportion of observations that were greater than predicted
quantiles = sapply(
  1:ncol(testing_observations),
  function(i){
    rowMeans(
      testing_observations[[i]] > unscale(t(future_observed[,,i])),
      na.rm = TRUE
    )
  }
)

# Observed versus nominal quantiles (should be on the 1:1 line)
plot(
  seq(0, 1, length = 101),
  sapply(seq(0, 1, length = 101),
         function(x) mean(quantiles < x, na.rm = TRUE)),
  type = "l",
  col = "purple",
  xlab = "nominal quantiles",
  ylab = "observed quantiles",
  xaxs = "i",
  yaxs = "i"
)
abline(0, 1, col = "#00000050")

# save output -------------------------------------------------------------
out = list(
  pt_est = data_frame(
    site_id = colnames(arima_data)[c(col(predicted_means))],
    model = "AR1SS",
    timeperiod = seq(last_training_year + 1, end_yr)[c(row(predicted_means))],
    obs = unlist(testing_observations),
    pt_fcast = c(predicted_means)
  ),
  intervals = data_frame(
    site_id = colnames(arima_data)[c(col(predicted_means))],
    model = "AR1SS",
    timeperiod = seq(last_training_year + 1, end_yr)[c(row(predicted_means))],
    obs = unlist(testing_observations),
    pt_fcast = c(predicted_means),
    lo = c(lower_CIs),
    hi = c(upper_CIs)
  ),
  quantiles = c(quantiles)
)
saveRDS(out, "stan-models/output.RDS")
