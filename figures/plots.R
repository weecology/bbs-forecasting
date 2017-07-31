devtools::load_all()
library(tidyverse)
library(cowplot)
library(ggjoy)
settings = yaml::yaml.load_file("settings.yaml")

timeframe = "train_22"
prepend_timeframe = function(x) {
  paste0("results", "/", timeframe, "/", x)
}


forecast_results = readRDS(prepend_timeframe("forecast.rds")) %>% 
  filter(!is.na(richness))
gbm_results = bind_rows(readRDS(prepend_timeframe("gbm_TRUE.rds")),
                        readRDS(prepend_timeframe("gbm_FALSE.rds")))
mistnet_results = lapply(dir(prepend_timeframe("mistnet_output/"), 
                             full.names = TRUE,
                             pattern = "TRUE|FALSE"),
                         readRDS) %>% 
  bind_rows()


average_results = bind_rows(readRDS(prepend_timeframe("avg_TRUE.rds")), 
                            readRDS(prepend_timeframe("avg_FALSE.rds")))

my_var = function(x){
  ifelse(length(x) == 1, 0, var(x))
}


rf_results = bind_rows(readRDS(prepend_timeframe("rf_predictions/all_FALSE.rds")),
                       readRDS(prepend_timeframe("rf_predictions/all_TRUE.rds")))

bound = bind_rows(forecast_results, gbm_results, average_results, 
                  rf_results, mistnet_results) %>% 
  group_by(site_id, year, richness, model, use_obs_model) %>% 
  summarize(sd = sqrt(mean(sd^2 + my_var(mean))), mean = mean(mean)) %>% 
  mutate(diff = richness - mean, z = diff / sd,
         p = pnorm(richness, mean, sd),
         deviance = -2 * dnorm(richness, mean, sd, log = TRUE)) %>% 
  ungroup()


# scatterplot -------------------------------------------------------------

ts_models = c("average", "naive", "auto.arima")
env_models = c("richness_gbm", "rf_sdm", "mistnet")
models = c(ts_models, env_models)

R2s = filter(bound, use_obs_model, year %in% c(2004, 2013)) %>% 
  mutate(x = min(mean), y = max(richness)) %>% 
  group_by(model, x, y, year) %>% 
  mutate(R2 = 1 - var(mean - richness) / var(richness)) %>% 
  mutate(R2_formatted = paste("R^2: ", 
                              format(R2,digits = 2, nsmall = 2)))

make_scatterplots = function(models){
  x_range = range(bound$mean)
  ggplot(data = filter(bound, use_obs_model, year %in% c(2004, 2013),
                       model %in% models), 
         aes(x = mean, y = richness)) +
    geom_hex() + 
    geom_abline(intercept = 0, slope = 1, alpha = .5) +
    facet_grid(year ~ model) +
    coord_equal() + 
    scale_fill_continuous(low = "gray90", high = "navy",
                          limits = c(1, 30)) +
    geom_text(inherit.aes = FALSE,
              data = filter(R2s, model %in% models),
              aes(x = x, y = y, label = R2_formatted),
              hjust = 0, vjust = 1,
              size = 3) +
    theme_light(base_size = 14) +
    xlim(x_range)
}

plot_grid(
  make_scatterplots(ts_models),
  make_scatterplots(env_models),
  align = "h",
  nrow = 2,
  labels = c("A. Time-series models", "B. Environmental models"))
ggsave(filename = "figures/scatter.png", width = 7.5, height = 10)

# Violins -----------------------------------------------------------------


for_violins = bound %>% 
  filter(model == "average", use_obs_model) %>% 
  left_join(filter(bound, model != "average", use_obs_model), 
            c("site_id", "year", "use_obs_model"),
            suffix = c(".avg", ".other")) %>% 
  mutate(dev_diff = deviance.other - deviance.avg,
         squared_diff = diff.other^2 - diff.avg^2,
         abs_diff = abs(diff.other) - abs(diff.avg))

make_violins = function(data, ylab = "", ylim, yintercept, adjust = 1, main){
  data %>% 
    ggplot(aes(x = model, y = y)) + 
    geom_hline(yintercept = yintercept, color = alpha(1, .25), size = 1/2) +
    geom_violin(fill = alpha("cornflowerblue", .9), size = 0, adjust = adjust) + 
    stat_summary(fun.data = "mean_cl_boot", colour = "darkblue", geom = "point",
                 size = 1) + 
    coord_cartesian(expand = FALSE, ylim = ylim) +
    theme_cowplot(9) +
    ylab(ylab) +
    xlab("Model type") +
    ggtitle(main)
}

# Note that the y axis doesn't extend all the way on this plot because
# of a small number of extreme outliers
plot_grid(
  for_violins %>% 
    mutate(model = model.other, y = -abs_diff) %>% 
    make_violins(yintercept = 0, adjust = 5,
                 ylim = quantile(for_violins$abs_diff, c(.0005, .9995)),
                 main = "A. Model improvement over\n\"Average\" baseline",
                 ylab = "Absolute error reduction"),
  for_violins %>% 
    mutate(model = model.other, y = plogis(-dev_diff / 2)) %>%
    make_violins(yintercept = 0.5, adjust = 1, ylim = c(0, 1), 
                 main = "B. Model posterior weight versus\n\"Average\" baseline",
                 ylab = "Posterior weight"),
  nrow = 2
)
ggsave(file = "figures/model_violins.png", width = 4, height = 5)


# Metrics over time -------------------------------------------------------

d = bound %>% 
  group_by(year, model, use_obs_model) %>% 
  summarize(rmse = sqrt(mean(diff^2)), 
            mean_deviance = -2 * mean(dnorm(richness, mean, sd, log = TRUE)),
            coverage = mean(p > .025 & p < .975)) %>% 
  filter(use_obs_model) %>% 
  gather(key = "variable", value = "value", rmse, mean_deviance, coverage) %>% 
  mutate(variable = forcats::fct_relevel(variable, "rmse", "mean_deviance", 
                                         "coverage"))

# Possibly-useful alternative 6-color colorblind palette
colors = c("#000000", dichromat::colorschemes$Categorical.12[c(4, 6, 10, 11, 12)])

ggplot() +
  geom_line(data = data_frame(x = c(min(d$year), max(d$year)), 
                              y = c(.95, .95),
                              variable = factor("coverage", 
                                                levels = levels(d$variable))),
            aes(x = x, y = y), inherit.aes = FALSE) +
  geom_line(data = d, aes(x = year, y = value, color = model), size = 1) + 
  scale_color_brewer(palette = "Set1") +
  ylab("") +
  facet_grid(variable~., scales = "free_y", switch = "y") +
  scale_x_continuous(expand = c(0, 0)) +
  theme_light(base_size = 10)
ggsave(file = "figures/performance_time.png", width = 3.75, height = 6)  


# digging into observer error ---------------------------------------------

obs_model = readRDS(prepend_timeframe("observer_model.rds"))

observer_uncertainties = obs_model$data %>% 
  filter(!in_train) %>% 
  group_by(observer_id, year, site_id, training_observer) %>% 
  summarize(observer_sd = sd(observer_effect)) %>% 
  ungroup() %>% 
  mutate(obs_class = cut_number(observer_sd, 10)) %>% 
  left_join(filter(bound, model == "average", use_obs_model), 
            by = c("year", "site_id"))

observer_uncertainties %>% 
  ggplot(aes(x = observer_sd^2, y = factor(year), 
             fill = year)) + 
  viridis::scale_fill_viridis(option = "B", guide = FALSE,
                              direction = 1) + 
  geom_joy(scale = 5, bandwidth = .75, color = "gray40") + 
  theme_joy(font_size = 11) +
  ylab(expression("Year" %->% "")) + 
  xlab("Observer-level uncertainty (squared error)") + 
  theme(axis.title.x = element_text(hjust = .5),
        axis.title.y = element_text(hjust = .5))
ggsave("figures/observer_uncertainty.png", width = 3.75, height = 4.5)


# Time series -------------------------------------------------------------


# sample_site_id = 88005
# sample_site_id = 91012
# sample_site_id = 2024
# sample_site_id = 91062
sample_site_id = 38032
sample_site_id = 72035

make_ts_plots = function(models, ylim, use_obs_model, sample_site_id, 
                         main = ""){
  grid = crossing(model = c(ts_models, env_models),
                  use_obs_model = unique(c(use_obs_model, FALSE)),
                  year = unique(obs_model$data$year))
  
  time_series_data = obs_model$data %>% 
    filter(site_id == sample_site_id, iteration == 1) %>% 
    select(site_id, year, richness, observer_id, richness) %>% 
    left_join(grid, by = "year") %>% 
    left_join(select(bound, -richness), by = c("year", "model", "use_obs_model",
                                               "site_id")) %>% 
    mutate(model = forcats::fct_relevel(model, unique(models)),
           observer_id = ifelse(use_obs_model, observer_id, 0)) %>% 
    select(site_id, year, model, use_obs_model, mean, sd, richness, observer_id,
           deviance)
  
  faceting = if (use_obs_model) {
    facet_grid(model ~ use_obs_model)
  } else {
    (facet_grid(~ model))
  }
  
  title = if (main != "") {
    ggtitle(main)
  } else {
    NULL
  }
  
  ggplot(filter(time_series_data, model %in% models), aes(x = year)) +
    geom_ribbon(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd),
                fill = "gray80") +
    geom_ribbon(aes(ymin = mean - 1 * sd, ymax = mean + 1 * sd),
                fill = "gray60") +
    geom_line(aes(y = mean), color = "gray30") +
    faceting +
    geom_line(aes(y = richness, alpha = .5 * (year < min(bound$year - 1)))) + 
    geom_point(aes(y = richness, fill = factor(observer_id)),
               size = 1, shape = 21, stroke = .5) +
    scale_alpha(range = c(0, 1), guide = FALSE) + 
    coord_cartesian(ylim = ylim, expand = TRUE) +
    ylab("Richness") +
    viridis::scale_fill_viridis(discrete = TRUE, option = "A", guide = FALSE) + 
    theme_light(base_size = 11) + 
    theme(plot.margin = unit(c(10, 7, 1, 7), units = "pt"),
          strip.background = element_rect(fill="white"), 
          strip.text = element_text(color = "black")) +
    title +
    xlab("Year")
}

# The warning about missing values is just saying that there are no predictions
# before 2004.
list(ts_models, env_models) %>% 
  map(~make_ts_plots(.x, 
                     ylim = c(42, 73), 
                     use_obs_model = FALSE, 
                     sample_site_id = 91034,
                     main = ifelse(
                       all(.x %in% ts_models),
                       "A. Time series models", 
                       "B. Environmental models"))
  ) %>% 
  plot_grid(plotlist = ., nrow = 2)
  
ggsave(filename = "figures/model_predictions.png", width = 7.5, height = 4)


make_ts_plots(c("average", "naive", "rf_sdm"), ylim = c(41, 91), 
              use_obs_model = TRUE, sample_site_id = 72035,
              main = "Controlling for observer differences")
ggsave(filename = "figures/observer_predictions.png", width = 4, height = 4)

# Numbers for the manuscript ----------------------------------------------


forecast_intercept_rate = forecast_results %>% 
  filter(model == "auto.arima", use_obs_model) %>% 
  pull(coef_names) %>% 
  map(~.x == "intercept") %>% 
  flatten_lgl() %>% 
  mean()
forecast_intercept_percent = round(100 * forecast_intercept_rate)


observer_deviance_data = bound %>% 
  mutate(use_obs_model = forcats::fct_recode(factor(use_obs_model), 
                                             no = "FALSE", 
                                             yes = "TRUE")) %>% 
  select(site_id, year, model, use_obs_model, deviance) %>% 
  spread(key = use_obs_model, value = deviance) %>% 
  mutate(y = plogis(0.5 * (no - yes)))

observer_error_data = bound %>% 
  mutate(use_obs_model = forcats::fct_recode(factor(use_obs_model), 
                                             no = "FALSE", 
                                             yes = "TRUE")) %>% 
  select(site_id, year, model, use_obs_model, diff) %>% 
  spread(key = use_obs_model, value = diff) %>% 
  mutate(y = abs(no) - abs(yes))


plot_grid(
  make_violins(observer_error_data, 
               main = "Absolute error reduction from observer model",
               ylab = "Reduction in absolute\nrichness error (species)",
               ylim = quantile(observer_error_data$y, c(.0005, .9995)), 
               yintercept = 0,
               adjust = 2),
  make_violins(observer_deviance_data, 
               main = "Posterior weight of model including observer effect",
               ylab = "Posterior weight of\nobserver model",
               ylim = c(0, 1), 
               yintercept = 0.5, 
               adjust = 5),
  nrow = 2
)
ggsave("figures/observers.png", width = 8, height = 5)


# Numerical odds and ends for manuscript ----------------------------------
bbs_data = get_pop_ts_env_data(
  start_yr = settings$start_yr, 
  end_yr = settings$timeframes[[timeframe]]$end_yr, 
  last_train_year = settings$timeframes[[timeframe]]$last_train_year, 
  min_year_percentage = settings$min_year_percentage)

variance_components = data_frame(site = obs_model$site_sigma^2, 
                                 observer = obs_model$observer_sigma^2,
                                 residual = obs_model$sigma^2)
variance_frac = colMeans(variance_components / rowSums(variance_components))

ms_numbers = list(
  N_species =  bbs_data %>% 
    distinct(species_id) %>% 
    nrow(),
  N_runs = obs_model$data %>% 
    distinct(site_id, year) %>% 
    nrow(),
  N_sites = obs_model$data %>% 
    distinct(site_id) %>% 
    nrow(),
  N_predict_sites = bound %>% 
    distinct(site_id) %>% 
    nrow(),
  N_predictions = bound %>% 
    distinct(site_id, year) %>% 
    nrow(),
  N_obs = obs_model$data %>% 
    distinct(observer_id) %>% 
    nrow(),
  richness_summary = obs_model$data %>% 
    distinct(site_id, year, richness) %>% 
    summarize(mean = round(mean(richness)),
              min = min(richness),
              max = max(richness)),
  n_gbm = gbm_results %>% 
    filter(site_id == 2001, year == 2004) %>% 
    group_by(use_obs_model) %>% 
    summarize(means = mean(n.trees)) %>% 
    pull(means) %>% 
    signif(2),
  rf_coverage_pct = bound %>% 
    filter(model == "rf_sdm", use_obs_model) %>% 
    summarize(round(100 * mean(p > .025 & p < .975))) %>% 
    pull(),
  all_R2s = R2s %>% 
    pull(R2) %>% 
    range() %>% 
    format(digits = 2) %>% 
    paste(collapse = "-"),
  env_R2s = R2s %>% 
    filter(model %in% env_models) %>% 
    pull(R2) %>% 
    range() %>% 
    format(digits = 2) %>% 
    paste(collapse = "-"),
  variance_percent = as.list(round(variance_frac * 100)),
  pct_is_average = forecast_results %>% 
    group_by(use_obs_model) %>% 
    filter(model == "auto.arima") %>% 
    summarize(is_avg = mean(100 * map_lgl(coef_names, ~all(.x == "intercept")))) %>% 
    pull(is_avg) %>% 
    round()
)  
cat(yaml::as.yaml(ms_numbers), file = "manuscript/numbers.yaml")

fitted = obs_model$data %>% 
  filter(in_train) %>% 
  group_by(site_id, year) %>% 
  summarize(diff = mean(richness) - 
              (mean(site_effect + observer_effect) + mean(obs_model$intercept))) %>% 
  ungroup()

ars = fitted %>% 
  complete(site_id, year) %>% 
  group_by(site_id) %>% 
  arrange(year) %>% 
  summarize(model = list(try(arima(diff, c(1, 0, 0),
                                   include.mean = FALSE))))
mean(ars$model %>% map_dbl(coef))



y = fitted %>% 
  complete(site_id, year) %>% 
  group_by(site_id) %>% 
  arrange(year) %>% spread(key = site_id, value = diff) %>% 
  dplyr::select(-year) %>% 
  as.matrix()

library(mvtnorm)

f = function(x) {
  b = x["b"]
  sigma = x["sigma"]
  covariance = matrix(NA, nrow = nrow(y), ncol = nrow(y))
  for (i in 1:nrow(y)) {
    for (j in 1:nrow(y)) {
      covariance[i, j] = sigma^2 * exp(-b * abs(i - j))
    }
  }
  - sum(
    sapply(
      1:ncol(y),
      function(i){
        non_na = !is.na(y[,i])
        dmvnorm(y[non_na,i], sigma = covariance[non_na, non_na], log = TRUE)
      }
    )
  )
}
o = optim(c(b = exp(-2), sigma = sd(y, na.rm = TRUE)), f, control = list(trace = 1))

exp(-o$par["b"])
