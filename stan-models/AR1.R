
if (as.numeric(packageDescription("forecast", fields = "Version")) < 7.1) {
  library("forecast", .libPaths()[[2]])
}
library(dplyr)
library(tidyr)
library(rstan)
devtools::load_all()

# MCMC chain count. Set to 4 for real analyses, set to 1 for quick ones
chains = 8

# Z-transform x based on the mean and sd of some other vector
rescale = function(x, baseline) {
  (x - mean(baseline)) / sd(baseline)
}

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25
last_training_year <- 2003
include_environment = 0


# The stan code expects the data to be sorted by site, then by year.
# Stan code currently can't handle NAs in environmental data, so we filter sites
# with missing data out.
if (!exists("raw", mode = "list")) {
  raw = get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
    group_by(site_id) %>%
    filter(length(lat) == end_yr - start_yr + 1) %>%
    ungroup() %>%
    complete(site_id, year) %>%
    arrange(site_id, year)
}

# Uncomment to drop most of the sites; useful for quick testing of the model
#raw = filter(raw, site_id %% 23 == 17)

# Split the data
training_observations = filter(raw, year <= last_training_year)
testing_observations = filter(raw, year > last_training_year)

# Build "data" object to pass to Stan -------------------------------------

# The autoregressive model in Stan will work on one long vector of richness
# estimates, one for each year/site combination (whether it was observed or
# not). Since the vector of observed richness values will be shorter than
# the full vector (because of missing values), we need a vector/index to
# map from one to the other.

data = training_observations

# give site/year combination a number in the observation index, then throw out
# the combinations that aren't used. This will be used to tell stan that ith
# row of the filtered matrix corresponds to the jth row of the original matrix.
data$observation_index = 1:nrow(data)
data = data[!is.na(data$richness), ]

# Convert observer_id to sequential integers to act as pointers in Stan.
data$observer_index = as.integer(factor(data$observer_id))

# Stan requires objects of different sizes in `data`, so rectangular tables
# won't work anymore
data = as.list(data)

# For all site/year combinations (including those with no observations),
# point to the corresponding site
data$N_train_years = length(unique(training_observations$year))
data$N_test_years = length(unique(testing_observations$year))
data$future_site_index = rep(seq_along(unique(data$site_id)),
                      each = data$N_test_years)

# Stan works better on z-scaled data.  Save the original mean and sd
# for later
data$scaled_richness = c(scale(data$richness))
center = mean(data$richness)
scale = sd(data$richness)

# Stan wants to know in advance how long the above vectors are going to be.
data$N_observations = length(data$scaled_richness)
data$N_sites = length(unique(data$site_id))
data$N_observers = max(data$observer_index)

# Stan won't accept data named "long" because it's a reserved word.
data$long = NULL

# We want *all* the predictor values, not just the years with data
data$env = scale(
  cbind(
    ndvi_sum = training_observations$ndvi_sum,
    log_elevs = log(training_observations$elevs)
  )
)
data$N_env = ncol(data$env)

data$future_env = cbind(
  ndvi_sum = rescale(testing_observations$ndvi_sum, training_observations$ndvi_sum),
  log_elevs = rescale(testing_observations$elevs, training_observations$elevs)
)

# If an observer has been seen before, use their observer_index. Otherwise,
# mark them as zero and pick a random value for them.
# There has to be a better way to do this. Maybe hadley/forcats to re-order
# factor levels before splitting into train and test?
data$future_observer_index = as.data.frame(data[c("observer_index", "observer_id")]) %>%
  distinct() %>%
  right_join(testing_observations, "observer_id") %>%
  arrange(site_id, year) %>%
  magrittr::extract2("observer_index")
data$future_observer_index[is.na(data$future_observer_index)] = 0

data$include_environment = include_environment

# Compile the stan model --------------------------------------------------

model = stan_model("stan-models/AR1.stan")

# Check the number of physical cores on the machine. With 2 or fewer cores,
# just use one core.  Otherwise, use all of them.
cores = ifelse(parallel::detectCores(logical = FALSE) <= 2,
               1,
               parallel::detectCores(logical = FALSE))


# Sample from the stan model ----------------------------------------------


samples = sampling(model, data = data, control = list(adapt_delta = 0.9),
                 cores = cores, chains = chains, refresh = 10, verbose = TRUE,
                 sample_file = "saved_sample", diagnostic_file = "diagnostic",
                 iter = 5000, thin = 2)

# Compress the stan output, delete the uncompressed versions
system("tar -czf tar -czf saved_sample*")
file.remove(dir("saved_sample[0-9]+.csv"))

# Save the model's predictions separately for easy loading
saveRDS(extract(samples, "future_observed")[[1]], "predictions.rds")

# extract model estimates -------------------------------------------------
unscale = function(data){
  data * scale + center
}
get_quantiles = function(x, q = c(.025, .975)){
  t(apply(x, 2, quantile, q))
}

extracted = rstan::extract(samples)

expected = structure(extracted$y, dim = c(nrow(extracted$y), data$N_train_years,
                                          data$N_sites))
future_expected = structure(extracted$future_y, dim = c(nrow(extracted$future_y), data$N_test_years,
                                          data$N_sites))

future_predicted = structure(extracted$future_observed, dim = c(nrow(extracted$future_observed), data$N_test_years,
                                                        data$N_sites))

# plot --------------------------------------------------------------------

k = sample.int(dim(expected)[3], 1)
matplot(start_yr:end_yr,
        rbind(unscale(t(expected[,,k])), unscale(t(future_expected[,,k]))),
        col = "#00000010", type = "l", lty = 1, ylab = "richness",
        xlab = "year", xlim = c(start_yr, end_yr),
        ylim = range(data$richness))
point_data = raw %>%
  filter(as.integer(factor(site_id)) == k) %>%
  mutate(observer_id = factor(observer_id))
point_data %>%
  lines(richness ~ year, data = ., col = 2, lwd = 2)
points(richness ~ year, data = point_data,
       col = 2 + as.integer(point_data$observer_id), pch = 16)
matlines(seq(last_training_year + 1, end_yr), t(apply(unscale(future_predicted[,,k]), 2, quantile, c(.025, .5, .975))), col = "purple", lty = 1, lwd = 2)
abline(v = last_training_year + 0.5, lty = 2)

# evaluation --------------------------------------------------------------
testing_observations$site_index = as.integer(factor(testing_observations$site_id))



f = function(g, ...){
  apply(unscale(extracted$future_observed), 2, g, ...)
}

results = data.frame(
  mean = f(mean),
  lower = f(quantile, .025),
  upper = f(quantile, .975),
  site_index = data$future_site_index,
  year = seq(last_training_year + 1, end_yr)
) %>%
  right_join(testing_observations, c("site_index", "year")) %>%
  mutate(diff = richness - mean, in_ci = richness > lower & richness < upper) %>%
  arrange(site_id, year) %>%
  mutate(quantile = rowMeans(t(unscale(extracted$future_observed)) > testing_observations$richness))




ggplot(NULL) +
  geom_violin(aes(x = results$year, y = abs(results$diff), group = factor(results$year))) +
  geom_smooth(aes(x = results$year, y = abs(results$diff)), se = FALSE)

plot(
  results %>% group_by(year) %>% summarize(coverage = mean(in_ci, na.rm = TRUE)),
  type = "o"
)
abline(h = 0.95)


# Observed versus nominal quantiles (should be on the 1:1 line)
plot(
  seq(0, 1, length = 101),
  sapply(seq(0, 1, length = 101),
         function(x) mean(results$quantile < x, na.rm = TRUE)),
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
    site_id = results$site_id,
    model = "AR1SS",
    timeperiod = results$year,
    obs = results$richness,
    pt_fcast = results$mean
  ),
  intervals = data_frame(
    site_id = results$site_id,
    model = "AR1SS",
    timeperiod = results$year,
    obs = results$richness,
    pt_fcast = results$mean,
    lo = results$lower,
    hi = results$upper
  ),
  quantiles = results$quantile
)
saveRDS(out, "stan-models/output.RDS")
