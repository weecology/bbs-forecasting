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

y1 = (x - mean(x)) / sd(x)
# Consider standardizing y1
data = list(N1 = nrow(y1), y1 = t(y1), N_sites = ncol(y1), N2 = 20)

m = stan_model(
  "gp/AR1.stan"
)

model = sampling(m, data = data, control = list(adapt_delta = 0.8))



# plot --------------------------------------------------------------------


extracted = rstan::extract(model)

k = sample.int(ncol(y1), 1)

N2 = 25
x = extracted$y[, k, nrow(y1)]
xx = matrix(x, 4000, N2)
mu = sigma = numeric(N2)
for(i in 1:N2){
  x = x * c(extracted$beta) + c(extracted$sigma) * rt(4000, extracted$nu)
  xx[,i] = x
  mu[i] = mean(x)
  sigma[i] = sd(x)
}

unscale = function(x){
  x * 9.5 + 55
}

add_nugget = function(x){
  out = rnorm(prod(length(x)), x, extracted$nugget)
  dim(out) = dim(x)
  out
}

N1 = nrow(y1)
N = N1 + N2


matplot(unscale(t(extracted$y[1:1000,k,])), col = "#00000020", xlim = c(1, N), ylim = range(arima_data, na.rm = TRUE), type = "n", lty = 1, ylab = "richness", xlab = "years")
points(unscale(y1[,k]), type = "o", pch = 16, lwd = 2, cex = 3/4, col = 2)
matlines(seq(1, N), t(apply(unscale(cbind(extracted$y[ ,k , ], xx)), 2, quantile, c(.025, .975))), col = 1, lty = 1, lwd = 2)
matlines(seq(N1 + 1, N), t(apply(unscale(add_nugget(xx)), 2, quantile, c(.025, .975))), col = 2, lty = 2, lwd = 2)
matlines(seq(N1 + 1, N), t(apply(unscale(add_nugget(xx)), 2, quantile, c(.005, .995))), col = 2, lty = 3, lwd = 2)
# legend("topleft", c("95% CI (\"true\")", "95% CI (observed)", "99% CI (observed)"), col = c(1, 2, 2),
#        lty = 1:3, lwd = 2)
