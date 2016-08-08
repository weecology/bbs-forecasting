library(forecast)
library(lme4)
library(rstan)
library(mvtnorm)
library(tidyr)
library(dplyr)
devtools::load_all()

start_yr <- 1982
end_yr <- 2013
min_num_yrs <- 25
last_training_year <- 2003
richness_w_env <- get_richness_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
  subset(year <=last_training_year) %>%
  na.omit()

arima_data = richness_w_env %>%
  dplyr::select(site_id, year, richness) %>%
  group_by(site_id) %>%
  tidyr::spread(key = year, value = richness) %>%
  t() %>%
  tbl_df()
colnames(arima_data) = arima_data[1, ]
arima_data = arima_data[-1, ]

x = as.matrix(arima_data[ , 0 == colSums(is.na(arima_data))])

indices = combn(1:20, 2)

z = t(
  sapply(1:ncol(indices),
         function(i){
           c(abs(indices[ , i][2] - indices[, i][1]), mean((x[indices[,i][1], ] - x[indices[,i][2], ])^2))
         }
  )
)
model = lm(z[,2] ~ z[,1])

plot(z, xlim = c(0, nrow(x)))
abline(model)
points(0, mean(apply(x, 2, var)), pch = 16)

summary(model)


# Stan --------------------------------------------------------------------

observed = (x - mean(x)) / sd(x)
data = list(N1 = nrow(observed), observed = t(observed), N_sites = ncol(observed), N2 = N2)

m = stan_model(
  "gp/AR1.stan"
)

model = sampling(m, data = data, control = list(adapt_delta = 0.8))



# plot --------------------------------------------------------------------


extracted = rstan::extract(model)

unscale = function(data){
  data * sd(x) + mean(x)
}

k = sample.int(ncol(x), 1)
matplot(
  seq(N1 + 1, N1+N2),
  unscale(t(apply(extracted$future_observed[,k,], 2, quantile, c(.005, .025, .975, .995)))),
  type = "l",
  lty = c(3, 2, 2, 3),
  col = 2,
  xlim = c(1, N1 + N2),
  ylim = range(x),
  ylab = "richness",
  xlab = "year",
  las = 1,
  lwd = 2
)
matlines(
  unscale(t(apply(cbind(extracted$y[,k,], extracted$future_y[,k,]), 2, quantile, c(.025, .975)))),
  col = "#00000070",
  lty = 1,
  lwd = 2
)
abline(v = N1 + 0.5, lty = 2, col = "#00000080", lwd = 2)
points(unscale(observed[, k]), cex = 3/4, pch = 16)
