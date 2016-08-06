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

stop()

final_values = extracted$y[ , , dim(extracted$y)[3]]
N2 = 25
CIs = lapply(
  1:N2,
  function(i){
    samples = structure(
      rnorm(length(final_values), mean = final_values, sd = sqrt(extracted$sigma^2 * i + extracted$nugget^2)),
      dim = dim(final_values)
    )
    out = cbind(
      tbl_df(t(apply(samples, 2, quantile, c(.025, .5, .975)))),
      site_num = 1:ncol(samples)
    )
    cbind(
      gather(out, key = "interval_type", value = value, -site_num),
      t_plus = i
    )
  }
) %>% bind_rows()

N1 = nrow(x)

i = sample.int(ncol(x), 1)
ggplot() +
  geom_point(data = data.frame(x = 1:N1, y = x[,i]), aes(x, y)) +
  geom_line(data = filter(CIs, site_num == i), aes(x = I(N1 + t_plus), y = value, group = interval_type)) +
  ylim(0, 100)
