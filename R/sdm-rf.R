#' @importFrom randomForest randomForest
one_rf_tree = function(bbs, vars, sp_id, iter, obs_model){
  d = bbs %>% 
    select(site_id, year, species_id, abundance) %>% 
    filter(species_id == sp_id) %>% 
    distinct() %>% 
    right_join(x_richness, by = c("site_id", "year"))
  
  my_formula = as.formula(paste("present ~", 
                                paste(vars, collapse = "+")))
   
  rf = randomForest(
    my_formula, 
    d %>% 
      filter(in_train, iteration == iter) %>% 
      mutate(present = factor(ifelse(is.na(abundance), 0, 1))),
    ntree = 1
  )  
  
  test = filter(d, !in_train, iteration == iter)
  test$mean = predict(rf, test, type = "prob")[,2]
  test$obs_model = obs_model
  
  select(test, site_id, year, species_id, mean, richness, obs_model, iteration)
}

rf_predict_species = function(sp_id, bbs, settings, x_richness, obs_model){
  vars = c(settings$vars, if (obs_model) {"observer_effect"})
  iters = sort(unique(x_richness$iteration))
  results = lapply(iters, one_rf_tree,
                   bbs = bbs, sp_id = sp_id, obs_model = obs_model, 
                   vars = vars) %>% 
    bind_rows()

  path = paste0("rf_predictions/sp_", sp_id, "_", obs_model, ".csv.gz")
  dir.create("rf_predictions", showWarnings = FALSE)
  write.csv(results, file = gzfile(path), row.names = FALSE)
  
  results %>% 
    group_by(site_id, year, species_id, richness) %>% 
    summarize(mean = mean(mean))
}

rf_predict_richness = function(bbs, x_richness, settings, obs_model) {
  
  parallel::mclapply(
    unique(bbs$species_id), 
    function(sp_id){
      rf_predict_species(sp_id, bbs = bbs, x_richness = x_richness, 
                         settings = settings, 
                         obs_model = obs_model)
    },
    mc.cores = 16,
    mc.preschedule = FALSE
  ) %>% 
    bind_rows() %>% 
    group_by(site_id, year, richness) %>% 
    summarize(sd = sqrt(sum(mean * (1 - mean))), mean = sum(mean)) %>% 
    ungroup() %>% 
    mutate(model = "rf_sdm", obs_model = obs_model)
}
