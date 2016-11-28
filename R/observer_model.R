#' Estimate observer effects
#' 
#' By default, \code{stan_file} points to a linear mixed model
fit_observer_model = function(stan_file = "mixed.stan", seed = 1,
                              output_file = "observer_model.rds"){
  library(rstan) # prevent annoying NAMESPACE issue described in
                 # https://github.com/stan-dev/rstan/issues/353
  
  set.seed(seed)
  settings = yaml::yaml.load_file("settings.yaml")
  
  d = get_pop_ts_env_data(settings$start_yr, settings$end_yr, 
                          settings$min_num_yrs) %>% 
    filter(!is.na(species_id))
  
  # Map arbitrary values of x to integers from 1:length(unique(x)),
  # used for indexing in the Stan model and in get_model_data below.
  make_index = function(x){
    as.integer(as.factor(x))
  }
  
  # Create a Stan-friendly version of the data.
  stan_df = d %>% 
    collapse_to_richness() %>% 
    filter(year <= settings$last_train_year) %>% 
    select(site_id, observer_id, richness, year) %>% 
    mutate(observer_index = make_index(observer_id), 
           site_index = make_index(site_id))
  
  # Append some information about the number of values to pass to stan
  stan_data = c(as.list(stan_df), 
                N = length(stan_df[[1]]), 
                N_observer = length(unique(stan_df$observer_id)), 
                N_site = length(unique(stan_df$site_id)))
  
  # Compile the model
  model = rstan::stan_model(stan_file)
  
  # Fit the model
  samples = rstan::sampling(model, chains = 1, data = stan_data)
  
  # Store the model output
  obs_model = c(
    rstan::extract(samples),
    stan_df,
    data = list(d)
  )
  
  
  saveRDS(obs_model, file = output_file)
  
  NULL
}

get_model_data = function(obs_model, sample_num){
  # One effect per site
  site_df = tibble::data_frame(
    site_id = unique(obs_model$site_id),
    site_effect = obs_model$site_effect[sample_num, unique(obs_model$site_index)]
  )
  
  # one effect per site-year combination, where data is available
  observer_df = tibble::data_frame(
    site_id = obs_model$site_id,
    observer_effect = obs_model$observer_effect[sample_num, obs_model$observer_index],
    expected_richness = obs_model$expected_richness[sample_num, ],
    year = obs_model$year
  )
  
  out = obs_model$data %>% 
    left_join(site_df, c("site_id")) %>% 
    left_join(observer_df, c("site_id", "year"))
  
  # Sanity check to ensure everything is in the right row:
  # Mean absolute difference between expected_richness + observer effect
  # and richness should be essentially zero
  mae = out %>% 
    collapse_to_richness() %>% 
    mutate(error = expected_richness + observer_effect - richness) %>% 
    summarize(mean(abs(error), na.rm = TRUE)) %>% 
    unlist()
  stopifnot(mae < 1E-9)
  
  out
}


