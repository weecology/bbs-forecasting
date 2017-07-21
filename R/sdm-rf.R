#' @importFrom randomForest randomForest
one_rf_tree = function(bbs, vars, sp_id, iter, use_obs_model, x_richness,
                       last_train_year, future, observer_sigmas,
                       settings){
  d = bbs %>% 
    select(site_id, year, species_id, abundance) %>% 
    filter(species_id == sp_id) %>% 
    distinct() %>% 
    right_join(filter(x_richness, iteration == iter), 
               by = c("site_id", "year")) %>% 
    mutate(present = factor(ifelse(is.na(abundance), 0, 1)))
  
  my_formula = as.formula(paste("present ~", 
                                paste(vars, collapse = " + ")))
  
  rf = randomForest(
    my_formula, 
    filter(d, in_train),
    ntree = 1
  )  
  
  test = make_test_set(d, future, observer_sigmas, settings) %>% 
    mutate(species_id = sp_id)
  test$mean = predict(rf, test, type = "prob")[,2]
  test$use_obs_model = use_obs_model
  
  select(test, site_id, year, species_id, mean, richness, use_obs_model, iteration)
}

rf_predict_species = function(sp_id, bbs, settings, x_richness, use_obs_model,
                              future, observer_sigmas, rf_dir){
  vars = c(settings$vars, if (use_obs_model) {"observer_effect"})
  iters = sort(unique(x_richness$iteration))
  results = lapply(iters, one_rf_tree,
                   bbs = bbs, sp_id = sp_id, use_obs_model = use_obs_model, 
                   vars = vars,
                   x_richness = x_richness,
                   last_train_year = settings$last_train_year,
                   future = future, 
                   observer_sigmas = observer_sigmas,
                   settings = settings) %>% 
    bind_rows()
  
  filename = paste0(rf_dir, "/sp_", sp_id, "_", use_obs_model, ".csv.gz")
  write.csv(results, file = gzfile(filename), row.names = FALSE)
  
  results %>% 
    group_by(site_id, year, species_id, richness, use_obs_model) %>% 
    summarize(mean = mean(mean))
}

rf_predict_richness = function(bbs, x_richness, settings, use_obs_model,
                               future, observer_sigmas, rf_dir) {
  
  out = parallel::mclapply(
    unique(bbs$species_id), 
    function(sp_id){
      rf_predict_species(sp_id, bbs = bbs, x_richness = x_richness, 
                         settings = settings, 
                         use_obs_model = use_obs_model,
                         future = future, 
                         observer_sigmas = observer_sigmas, rf_dir)
    },
    mc.cores = 8,
    mc.preschedule = FALSE
  ) %>% 
    purrr::map(combine_sdm_iterations) %>% 
    bind_rows() %>% 
    group_by(site_id, year, richness, use_obs_model) %>% 
    summarize(mean = sum(mean), sd = sqrt(sum(sd^2))) %>% 
    ungroup() %>% 
    mutate(model = "rf_sdm")
  
  out
}


combine_sdm_iterations = function(d){
  d %>% 
    group_by(site_id, year, species_id, richness, use_obs_model) %>% 
    summarize(mean = mean(mean), sd = sqrt(mean(mean * (1 - mean))))
}
