library(dplyr)
library(tidyr)
library(rstan)
devtools::load_all()

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25
last_training_year <- 2003

# The stan code expects the data to be sorted by site, then by year
raw = get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
  complete(site_id, year) %>%
  arrange(site_id, year)

training_observations = filter(raw, year <= last_training_year)
testing_observations = filter(raw, year > last_training_year)

n_train_years = length(unique(training_observations$year))
n_test_years = length(unique(testing_observations$year))


# Build "data" object to pass to Stan -------------------------------------

# The autoregressive model in Stan will work on one long vector of richness
# estimates, one for each year/site combination (whether it was observed or
# not). In the autoregressive model, each element that's not from year
# `start_yr` will point back to the previous element. For each elment, we also
# need to point to the site where the observation was taken and (if richness
# data was collected) to the observer_id. Since the vector of observed richness
# values will be shorter than the full vector (because of missing values), we
# also need a vector/index to map from one to the other.

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
data$N_train_years = length(start_yr:last_training_year)
data$site_index = rep(seq_along(unique(data$site_id)),
                      each = data$N_train_years)

# Point to all the site/year combinations (including those with no observations)
# that occurred in the first year of data collection, and separately to those
# that occurred after the first year of data collection.
year_mod = seq_along(data$site_index) %% data$N_train_years
data$which_first = which(year_mod == 1)
data$which_non_first = which(year_mod != 1)

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

# Compile the stan model --------------------------------------------------

m = stan_model("stan-models/AR1.stan")
stop()
# Check the number of physical cores on the machine. With 2 or fewer cores,
# just use one core.  Otherwise, use all of them.
cores = ifelse(parallel::detectCores(logical = FALSE) <= 2, 1,
               parallel::detectCores(logical = FALSE))


# Sample from the stan model ----------------------------------------------

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

expected = structure(extracted$y, dim = c(nrow(extracted$y), data$N_train_years,
                                          data$N_sites))

# Predict -----------------------------------------------------------------

yy = expected[ , dim(expected)[2], ]
future_y = array(NA, c(nrow(yy), ncol(yy), n_test_years))
for (i in 1:n_test_years) {
  future_means = extracted$site_alphas + extracted$site_betas * yy
  yy = rnorm(length(yy), future_means, extracted$sigma_autoreg)
  future_y[ , , i] = yy
}

# plot --------------------------------------------------------------------

k = sample.int(dim(expected)[3], 1)
predicted = unscale(t(future_y[,k,]))
matplot(start_yr:end_yr,
        rbind(unscale(t(expected[,,k])), predicted),
        col = "#00000010", type = "l", lty = 1, ylab = "richness",
        xlab = "year", xlim = c(start_yr, end_yr),
        ylim = range(unscale(data$observed)))
point_data = training_observations %>%
  filter(as.integer(factor(site_id)) == k) %>%
  mutate(observer_id = factor(observer_id))
point_data %>%
  lines(richness ~ year, data = ., col = 2, lwd = 2)
point_data %>%
  points(richness ~ year, data = ., col = 2 + as.integer(observer_id),
         pch = 16, ylim = range(unscale(data$observed)))

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
