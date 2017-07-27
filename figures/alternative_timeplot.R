i = i + 1
sample_site_id = rev(unique(bound$site_id))[i]
#sample_site_id = 14198
#sample_site_id = 92017
#sample_site_id = 91062
#sample_site_id = 91034
sample_site_id = 91022

grid = crossing(model = c(ts_models, env_models),
                use_obs_model = c(TRUE, FALSE),
                year = unique(obs_model$data$year))

time_series_data = obs_model$data %>% 
  filter(site_id == sample_site_id, iteration == 1) %>% 
  select(site_id, year, richness, observer_id, richness) %>% 
  left_join(grid, by = "year") %>% 
  left_join(select(bound, -richness), by = c("year", "model", "use_obs_model",
                                             "site_id")) %>% 
  mutate(model = forcats::fct_relevel(model, unique(models)),
         observer_id = ifelse(use_obs_model, observer_id, 0)) %>% 
  select(site_id, year, model, use_obs_model, mean, sd, richness, observer_id)


observer_colors = RColorBrewer::brewer.pal(4, "PuOr")

# time_series_data$observer_id = factor(time_series_data$observer_id)
# levels(time_series_data$observer_id) = c("0", "A", "B", "C")
ggplot(time_series_data, aes(x = year)) +
  geom_ribbon(aes(ymin = mean - 1.96 * sd, ymax = mean + 1.96 * sd,
                  group = use_obs_model,
                  alpha = .2 * as.numeric(!use_obs_model)),
              fill = "black") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd,
                  group = use_obs_model, 
                  alpha = .2 * as.numeric(!use_obs_model)),
              fill = "black") +
  geom_line(aes(y = mean, group = use_obs_model, color = use_obs_model,
                alpha = .5 * as.numeric(!use_obs_model))) +
  facet_wrap(~ model) +
  geom_line(aes(y = richness, alpha = .5 * (year < min(bound$year - 1)))) + 
  geom_point(aes(y = richness), size = 1, shape = 21, fill = "red") +
  scale_alpha(range = c(0, 1), guide = FALSE) + 
  coord_cartesian(xlim = c(1980, 2014), ylim = c(44, 83), expand = FALSE) +
  ylab("Richness") +
  scale_fill_brewer(guide = FALSE) +
  scale_color_manual(values = c("black", observer_colors[1]), 
                     guide = FALSE) +
  theme_light(base_size = 11) + 
  theme(plot.margin = unit(c(12, 7, 1, 7), units = "pt"),
        panel.spacing = unit(1, "lines"),
        strip.background = element_rect(fill="white"), 
        strip.text = element_text(color = "black"))
ggsave("figures/alternative_predictions.png", width = 7.5, height = 4)
