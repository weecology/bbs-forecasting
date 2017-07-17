devtools::load_all()
library(tidyverse)
library(cowplot)
settings = yaml::yaml.load_file("settings.yaml")

timeframe = "train_22"
prepend_timeframe = function(x) {
  paste0("results", "/", timeframe, "/", x)
}


forecast_results = readRDS(prepend_timeframe("forecast.rds")) %>% 
  filter(!is.na(richness))
gbm_results = bind_rows(readRDS("gbm_TRUE.rds"),
                        readRDS("gbm_FALSE.rds"))
mistnet_results = lapply(dir(prepend_timeframe("mistnet_output/"), 
                             full.names = TRUE,
                             pattern = "TRUE|FALSE"),
                         readRDS) 


mean(map_dbl(mistnet_results, ~mean(.x$mean - .x$richness)))


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
         deviance = -2 * dnorm(richness, mean, sd, log = TRUE))

bound %>% 
  filter(model == "average") %>% 
  left_join(filter(bound, model != "average"), 
            c("site_id", "year", "use_obs_model")) %>% 
  mutate(y = deviance.y - deviance.x) %>% 
  filter(use_obs_model) %>% 
  ggplot(aes(x = model.y, y = y)) +
  ggforce::geom_sina(binwidth = .1, size = .25) +
  geom_hline(yintercept = 0)


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
  filter(use_obs_model) %>% 
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

obs_model = readRDS(prepend_timeframe("observer_model.rds"))

observer_uncertainties = obs_model$data %>% 
  filter(!in_train) %>% 
  group_by(observer_id, year, site_id) %>% 
  summarize(observer_sd = sd(observer_effect)) 

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
lmer_data = filter(bound, model == "average", use_obs_model) %>% 
  left_join(observer_uncertainties, by = c("site_id", "year")) %>% 
  mutate(squared_diff = diff^2)

# Variance associated with observer differences versus year
lmer(squared_diff ~ poly(observer_sd, 2) + (1|site_id) + year, 
     data = lmer_data) %>% 
  anova()
