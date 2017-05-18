devtools::load_all()
library(tidyverse)
library(cowplot)
settings = yaml::yaml.load_file("settings.yaml")

forecast_results = readRDS("forecast.rds") %>% 
  filter(!is.na(richness))
gbm_results = bind_rows(readRDS("gbm_TRUE.rds"),
                        readRDS("gbm_FALSE.rds"))
mistnet_results = lapply(dir("mistnet_output/", full.names = TRUE,
                             pattern = "TRUE|FALSE"),
                         readRDS) 


mean(map_dbl(mistnet_results, ~mean(.x$mean - .x$richness)))


average_results = bind_rows(readRDS("avg_TRUE.rds"), readRDS("avg_FALSE.rds"))

my_var = function(x){
  ifelse(length(x) == 1, 0, var(x))
}


rf_results = bind_rows(readRDS("rf_predictions/all_FALSE.rds"),
                       readRDS("rf_predictions/all_TRUE.rds"))

bound = bind_rows(forecast_results, gbm_results, average_results, 
                  rf_results, mistnet_results) %>% 
  group_by(site_id, year, richness, model, use_obs_model) %>% 
  summarize(sd = sqrt(mean(sd^2 + my_var(mean))), mean = mean(mean)) %>% 
  mutate(diff = richness - mean, z = diff / sd,
         p = pnorm(richness, mean, sd),
         deviance = -2 * dnorm(richness, mean, sd, log = TRUE))


for_violins = bound %>% 
  group_by(site_id, year, use_obs_model, model) %>% 
  summarize(mse = Metrics::mse(richness, mean), mean_deviance = mean(deviance)) %>% 
  arrange(use_obs_model) %>% 
  group_by(site_id, year, model) %>% 
  summarize(n = n(), diff_mse = mse[2] - mse[1], diff_deviance = mean_deviance[2] - mean_deviance[1])

for_violins %>% 
  ggplot(aes(x = model, y = diff_mse)) + 
  geom_violin(fill = 2, size = 0, alpha = .5, adjust = 2) +
  geom_hline(yintercept = 0, color = alpha(1, .25)) +
  stat_summary(fun.data = "mean_cl_boot", colour = "darkred", geom = "point") +
  coord_cartesian(ylim = quantile(for_violins$diff_mse, c(.025, .975)), 
                  expand = FALSE)


d = bound %>% 
  group_by(year, model, use_obs_model) %>% 
  summarize(mse = mean(diff^2), 
            mean_deviance = -2 * mean(dnorm(richness, mean, sd, log = TRUE)),
            coverage = mean(p > .025 & p < .975))

# Possibly-useful alternative 6-color colorblind palette
# c("#000000", dichromat::colorschemes$Categorical.12[c(4, 6, 10, 11, 12)])

d %>% 
  filter(!use_obs_model) %>% 
  gather(key = "variable", value = "value", mse, mean_deviance, coverage) %>% 
  mutate(variable = forcats::fct_relevel(variable, "mse", "mean_deviance", 
                                         "coverage")) %>% 
  ggplot(aes(x = year, y = value, color = model)) +
  geom_line(size = 1) + 
  scale_color_brewer(palette = "Set1") +
  facet_grid(variable~., scales = "free_y", switch = "y") +
  theme(strip.background = element_blank()) +
  ylab("")

plots = map(c("mse", "mean_deviance", "coverage"),
            ~(d %>% 
                filter(!use_obs_model) %>% 
                ggplot(aes_string(x = "year", y = .x, color = "model")) +
                geom_line(size = 1) + 
                scale_color_brewer(palette = "Set1") +
                ggtitle(.x) +
                coord_cartesian(ylim = if(.x == "coverage"){c(min(d$coverage, na.rm = TRUE) - .01, 1.0)}, 
                                expand = .x != "coverage") +
                if(.x=="coverage"){geom_hline(yintercept = .95, color = alpha("black", .25))}))

plot_grid(plotlist = plots, nrow = 3)



# digging into observer error ---------------------------------------------

obs_model = readRDS("observer_model.rds")
av_point = average_results %>% 
  group_by(year, site_id, richness) %>% 
  summarize(mean = mean(mean))

obs_model$data %>% 
  filter(!in_train) %>% 
  group_by(year, training_observer, observer_id) %>% 
  summarize(resid_var = var(observer_effect)) %>% 
  group_by(year, training_observer) %>% 
  summarize(mean_var = mean(resid_var)) %>% 
  ggplot(aes(year, mean_var, color = training_observer)) +
  geom_line()


obs_model$data %>% 
  filter(!in_train) %>% 
  group_by(year) %>% 
  summarize(`% new` = mean(!training_observer)) %>% 
  plot()

observer_uncertainties = obs_model$data %>% 
  filter(!in_train) %>% 
  group_by(observer_id, year, site_id, training_observer) %>% 
  summarize(observer_sd = sd(observer_effect)) %>% 
  right_join(av_point, by = c("year", "site_id")) %>% 
  mutate(squared_error = (richness - mean)^2, log = TRUE) %>% 
  ungroup()


plot(gam(log(squared_error) ~ s(year) + s(observer_sd), 
         data = observer_uncertainties))

breaks = hist(observer_uncertainties$observer_sd, breaks = 15, 
              plot = FALSE)$breaks
binned = observer_uncertainties %>% 
  mutate(`observer sd` = cut(observer_sd, breaks = breaks)) %>% 
  group_by(`observer sd`, year) %>% 
  summarize(mse = mean(squared_error), count = n()) %>% 
  group_by(year) %>% 
  mutate(proportion = count / sum(count))

binned %>% 
  ggplot(aes(year, `observer sd`, fill = mse)) + 
  geom_raster() +
  viridis::scale_fill_viridis() + 
  coord_cartesian(expand = FALSE) + 
  ggtitle("error by observer sd estimate & time")

binned %>% 
  ggplot(aes(year, `observer sd`, fill = proportion)) + 
  geom_raster() +
  viridis::scale_fill_viridis() + 
  coord_cartesian(expand = FALSE) + 
  ggtitle("observer sd estimates over time")

library(gbm_results)
g = gbm_results(squared_error ~ ordered(year) + observer_sd, 
                data = observer_uncertainties,
                distribution = "gaussian",
                var.monotone = c(1, 1),
                n.trees = 1E4,
                shrinkage = 5E-4,
                interaction.depth = 3)
n.trees = gbm.perf(g)
plot(g, c("ordered(year)", "observer_sd"), 
     n.trees = n.trees,
     continuous.resolution = 100,
     return.grid = TRUE) %>% 
  rename(mse = y) %>% 
  ggplot(aes(`ordered(year)`, observer_sd, fill = mse)) +
  geom_raster() +
  viridis::scale_fill_viridis(option = "A", trans = "log",
                              breaks = c(25, 50, 75)) +
  coord_cartesian(expand = FALSE)


summary(g, n.trees = n.trees)
