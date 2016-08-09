library(dplyr)
library(tidyr)
library(rstan)
devtools::load_all()

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25
last_training_year <- 2003
richness_w_env <- get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
  na.omit()

arima_data = richness_w_env %>%
  dplyr::select(site_id, year, richness) %>%
  group_by(site_id) %>%
  spread(key = year, value = richness) %>%
  t() %>%
  tbl_df()

# First row is actually column names
colnames(arima_data) = arima_data[1, ]
arima_data = arima_data[-1, ]

training_observations = arima_data[1:length(start_yr:last_training_year), ]
testing_observations = arima_data[-(1:length(start_yr:last_training_year)), ]

# one long vector of non-NA observations
observed = as.matrix(arima_data)[!is.na(arima_data)]

# sites associated with each value of `observed`
site_id = c(col(arima_data)[!is.na(arima_data)])

# Number of observations assocaited with each site
site_length = rle(site_id)$lengths

# Which values are not missing?
index = seq(1, prod(dim(arima_data)))[!is.na(arima_data)]

# Stan --------------------------------------------------------------------

data = list(
  N_obs = length(observed),
  N1 = length(start_yr:last_training_year),
  N2 = length(seq(last_training_year + 1, end_yr)),
  N_sites = ncol(arima_data),
  index = index,
  observed = c(scale(observed))
)

m = stan_model(
  "gp/AR1.stan"
)

model = sampling(m, data = data, control = list(adapt_delta = 0.8))



# extract model estimates -------------------------------------------------
unscale = function(data){
  data * sd(observed) + mean(observed)
}
get_quantiles = function(x, q = c(.025, .975)){
  t(apply(x, 2, quantile, q))
}

extracted = rstan::extract(model)

y = array(extracted$y, dim = c(nrow(extracted$y), dim(arima_data)))
future_y = array(extracted$future_y, dim = c(nrow(extracted$future_y), data$N2, data$N_sites))
future_observed = array(extracted$future_observed, dim = c(nrow(extracted$future_y), data$N2, data$N_sites))

predicted_means = unscale(apply(future_y, 2:3, mean))
lower_CIs = unscale(apply(future_observed, 3, get_quantiles, .025))
upper_CIs = unscale(apply(future_observed, 3, get_quantiles, .975))

# plot --------------------------------------------------------------------

plot_site = function(k){
  matplot(
    seq(data$N1 + 1, data$N1 + data$N2),
    unscale(get_quantiles(future_observed[,,k], c(.005, .025, .975, .995))),
    type = "l",
    lty = c(3, 2, 2, 3),
    col = 2,
    xlim = c(1, data$N1 + data$N2),
    xlab = "year",
    ylab = "species richness",
    las = 1,
    lwd = 2,
    ylim = range(observed)
  )
  matlines(
    unscale(get_quantiles(cbind(y[,,k], future_y[,,k]))),
    lty = 1,
    lwd = 2,
    col = "#00000060"
  )
  points(arima_data[ , k], pch = 16, cex = 0.75)
  abline(v = data$N1 + 1/2, col = "#00000060", lty = 2, lwd = 2)
  lines(unscale(colMeans(cbind(y[,,k], future_y[,,k]))))
}

par(mfrow = c(2, 2))
sapply(sample.int(ncol(arima_data), 4), plot_site)
par(mfrow = c(1, 1))


# evaluation --------------------------------------------------------------

diffs = as.matrix(predicted_means - testing_observations)

# RMSE
hist(diffs)

# Proportion outside the confidence interval (should each be close to 0.025)
mean(testing_observations < lower_CIs, na.rm = TRUE) # below lower
mean(testing_observations > upper_CIs, na.rm = TRUE) # above upper
