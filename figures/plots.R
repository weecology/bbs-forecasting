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

R2s = filter(bound, use_obs_model, year %in% c(2004, 2013)) %>% 
  mutate(x = min(mean), y = max(richness)) %>% 
  group_by(model, x, y, year) %>% 
  summarize(R2 = paste("R^2: ", 
                       format(1 - var(mean - richness) / var(richness),
                              digits = 2, nsmall = 2)))
            
ggplot(data = filter(bound, use_obs_model, year %in% c(2004, 2013)), 
       aes(x = mean, y = richness)) +
  geom_hex() + 
  geom_abline(intercept = 0, slope = 1, alpha = .5) +
  facet_grid(year ~ model) +
  coord_equal() + 
  scale_fill_continuous(low = "gray90", high = "navy") +
  geom_text(inherit.aes = FALSE,
            data = R2s,
            aes(x = x, y = y, label = R2),
            hjust = 0, vjust = 1,
            size = 3) +
  theme_light()


# Violins -----------------------------------------------------------------


for_violins = bound %>% 
  filter(model == "average", use_obs_model) %>% 
  left_join(filter(bound, model != "average"), 
            c("site_id", "year", "use_obs_model"),
            suffix = c(".avg", ".other")) %>% 
  mutate(dev_diff = deviance.other - deviance.avg,
         abs_diff = abs(diff.other) - abs(diff.avg))

for_violins %>% 
  ggplot(aes(x = model.other, y = -abs_diff)) + 
  geom_violin(fill = 2, size = 0, alpha = .5, adjust = 1) +
  geom_hline(yintercept = 0, color = alpha(1, .25)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "darkred", geom = "point") +
  coord_cartesian(expand = FALSE) +
  ylab("Improvement over 'Average' (absolute error)")

for_violins %>% 
  ggplot(aes(x = model.other, y = plogis(dev_diff / 2))) + 
  geom_violin(fill = 2, size = 0, alpha = .5, adjust = .5) +
  geom_hline(yintercept = 0.5, color = alpha(1, .25)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "darkred", geom = "point") +
  coord_cartesian(expand = FALSE) +
  ylim(0, 1)



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
# c("#000000", dichromat::colorschemes$Categorical.12[c(4, 6, 10, 11, 12)])

theme_gray_elements =  theme(panel.background = theme_gray()$panel.background,
                             panel.grid.major = theme_gray()$panel.grid.major,
                             panel.grid.minor = theme_gray()$panel.grid.minor)

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
  theme_gray_elements
  


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
  ggplot(aes(x = year, y = sqrt(diff^2), group = obs_class, 
             color = as.integer(obs_class))) + 
  geom_smooth(alpha = 0, method = "lm") + 
  viridis::scale_color_viridis() +
  geom_smooth(aes(group = NULL), size = 3, color = "black")

observer_uncertainties %>% 
  ggplot(aes(x = observer_sd^2, y = factor(year), 
             fill = year)) + 
  viridis::scale_fill_viridis(option = "B", guide = FALSE,
                              direction = 1) + 
  geom_joy(scale = 5, bandwidth = .75, color = "gray40") + 
  theme_joy() +
  ylab(expression("Year" %->% "")) + 
  xlab("Observer-level uncertainty") + 
  theme(axis.title.x = element_text(hjust = .5),
        axis.title.y = element_text(hjust = .5))

library(lme4)
# Variance associated with observer differences versus year
lmer(I(diff^2) ~ poly(observer_sd, 2) + (1|site_id) + year, 
     data = observer_uncertainties) %>% 
  anova()

library(mgcv)
bs = smooth.construct(s(observer_sd), observer_uncertainties, NULL)$X
l = lmer(scale(diff^2) ~ bs[,-9] + (1|site_id) + scale(year), 
     data = observer_uncertainties)

observer_uncertainties %>% 
  group_by(year) %>% 
  summarize(mse = mean((richness - mean)^2), moe = mean(observer_sd^2)) %>% 
  gather(key = var, value = y, mse, moe) %>% 
  ggplot(aes(x = year, y = y, color = var)) + 
  geom_line()


# Time series -------------------------------------------------------------

#sample_site_id = unique(bound$site_id)[i]
#sample_site_id = 82006
#sample_site_id = 88005
#sample_site_id = 2022
#sample_site_id = 91012
#sample_site_id = 2024
sample_site_id = 82003
sample_site_id = sample(unique(bound$site_id), 1)

models = c("average", "naive", "auto.arima")
grid = crossing(model = models,
                use_obs_model = c(TRUE, FALSE),
                year = unique(obs_model$data$year))

time_series_data = obs_model$data %>% 
  filter(site_id == sample_site_id, iteration == 1) %>% 
  select(site_id, year, richness, observer_id, richness) %>% 
  left_join(grid, by = "year") %>% 
  left_join(select(bound, -richness), by = c("year", "model", "use_obs_model",
                                             "site_id")) %>% 
  mutate(model = forcats::fct_relevel(model, models),
         observer_id = ifelse(use_obs_model, observer_id, 0)) %>% 
  select(site_id, year, model, use_obs_model, mean, sd, richness, observer_id)


observer_colors = RColorBrewer::brewer.pal(4, "PuOr")

# The warning about missing values is just saying that there are no predictions
# before 2004.
ggplot(time_series_data, aes(x = year)) +
  geom_ribbon(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd),
              fill = "gray80") +
  geom_ribbon(aes(ymin = mean - 1 * sd, ymax = mean + 1 * sd),
              fill = "gray60") +
  geom_line(aes(y = mean), size = 1.5, color = "gray30") +
  facet_grid(use_obs_model ~ model) +
  geom_line(aes(y = richness, alpha = as.numeric(year < min(bound$year - 1)))) + 
  geom_point(aes(y = richness, fill = factor(observer_id)),
             size = 3.5, shape = 21) +
  scale_alpha(range = c(0, 1), guide = FALSE) + 
  coord_cartesian(ylim = c(35, 70), expand = TRUE) +
  ylab("Richness") +
  scale_fill_manual(values = c("black", observer_colors), guide = FALSE) +
  theme(panel.grid.major = theme_light()$panel.grid.major,
        panel.grid.minor = theme_light()$panel.grid.minor)

forecast_results %>% filter(site_id == sample_site_id, model == "auto.arima", use_obs_model) %>% pull(coef_names)


# Auto.arima --------------------------------------------------------------

forecast_results %>% 
  filter(model == "auto.arima", !use_obs_model) %>% 
  pull(coef_names) %>% 
  map_chr(paste, collapse = " ") %>% 
  table() %>% 
  sort()

forecast_results %>% 
  filter(model == "auto.arima", use_obs_model) %>% 
  pull(coef_names) %>% 
  map_chr(paste, collapse = " ") %>% 
  table() %>% 
  sort()
