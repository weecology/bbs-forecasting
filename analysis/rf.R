library(rpart)
library(tidyverse)
devtools::load_all()

settings = yaml::yaml.load_file("settings.yaml")
iters = 1:200 # During interactive use/testing, only include a few MCMC samples


obs_model = readRDS("observer_model.rds")
x_richness = obs_model$data %>% 
  mutate(intercept = c(obs_model$intercept[iteration]), 
         sd = c(obs_model$sigma[iteration]),
         expected_richness = richness - observer_effect) %>% 
  filter(iteration %in% iters)

bbs = get_bbs_data() %>% 
  select(site_id, year, species_id) %>% 
  mutate(presence = 1)


sdm_tree = function(spid, iter, obs_model){
  d = x_richness %>% 
    filter(iteration == iter) %>% 
    left_join(filter(bbs, species_id == spid), c("year", "site_id")) %>% 
    mutate(presence = ifelse(is.na(presence), 0, 1))
  
  vars = c(settings$vars, "presence", if (obs_model) {"observer_effect"})
  
  boot_rows = sample.int(sum(d$in_train), replace = TRUE)
  non_boot_rows = seq(1, sum(d$in_train)) %>% discard(~.x %in% boot_rows)
  
  train = d %>% filter(in_train) %>% .[boot_rows, vars] # Bootstrap resample
  oob = d %>% filter(in_train) %>% .[non_boot_rows, vars] # Out of bag
  test = d %>% filter(!in_train)
  
  if (var(train$presence) == 0) {
    # if presence is constant, return that constant
    p = rep(train$presence[1], nrow(test))
  } else{
    t = rpart(presence ~ ., data = train[boot_rows, ], method = "class",
              maxcompete = 0, maxsurrogate = 0, cp = 0.00001, minbucket = 25,
              xval = 0)
    p = predict(t, test)[,2]
  }
  
  data_frame(
    iteration = iter,
    species_id = spid,
    site_id = test$site_id,
    year = test$year,
    obs_model = obs_model,
    p = p
  )
}


sdm = function(iters, spid, obs_model) {
  map(iters, ~sdm_tree(spid = spid, iter = .x, obs_model = TRUE)) %>% 
    bind_rows() %>% 
    group_by(species_id, site_id, year, obs_model) %>% 
    mutate(variance = p * (1 - p)) %>% 
    summarize(mean = mean(p), sd = sqrt(mean(variance) + var(p))) %>% 
    ungroup()
}


sdms = map(
  unique(bbs$species_id),
  ~sdm(iters = iters, .x, TRUE)
) %>% 
  bind_rows()

sdms_no_obs = map(
  unique(bbs$species_id),
  ~sdm(iters = iters, .x, FALSE)
) %>% 
  bind_rows()



richness_rf = function(d, obs_model, settings){
  if (obs_model) {
    d$y = d$expected_richness
  } else {
    d$y = d$richness
  }
  
  boot_rows = sample.int(sum(d$in_train), replace = TRUE)
  non_boot_rows = seq(1, sum(d$in_train)) %>% discard(~.x %in% boot_rows)
  
  train = d %>% filter(in_train) %>% .[boot_rows, ] # Bootstrap resample
  oob = d %>% filter(in_train) %>% .[non_boot_rows, ] # Out of bag
  test = d %>% filter(!in_train)
  
  # First column is the intercept; drop it
  x = cbind(
    model.matrix(as.formula(settings$formula), train)[, -1],
    y = train$y
  )
  
  t = rpart(y ~ ., data = as.data.frame(x), method = "anova", 
            maxcompete = 0, maxsurrogate = 0, cp = 0.00001, minbucket = 25,
            xval = 0)
  
  # Store the SD associated with each leaf of the tree (i.e. each `mean` value)
  sd_table = data_frame(y = x[,"y"], mean = predict(t)) %>% 
    group_by(mean) %>% 
    summarize(sd = sqrt(mean((y - mean)^2)))
  
  test = test %>% 
    mutate(mean = predict(t, test)) %>% 
    select(-sd) %>% 
    left_join(sd_table, "mean")
  
  if (obs_model) {
    test$mean = test$mean + test$observer_effect
  }
  
  test %>% 
    mutate(model = "rf", obs_model = obs_model) %>% 
    select(site_id, year, mean, sd, richness, model, obs_model)
}

rf_predictions = x_richness %>% 
  group_by(iteration) %>% 
  by_slice(richness_rf, obs_model = TRUE, settings = settings, .collate = "rows")

rf_no_obs = x_richness %>% 
  group_by(iteration) %>% 
  by_slice(richness_rf, obs_model = FALSE, settings = settings, .collate = "rows")

