#' Estimate observer effects
#' 
#' The default `mixed.stan` file is a linear mixed model for richness
#' using site_id and observer_id as random effects.
#' 
#' @return \code{NULL}. Saves a list of relevant values in `output_file`
#'   and the `stanfit` object to disk.
fit_observer_model = function(stan_file = "mixed.stan", seed = 1,
                              output_file = "observer_model.rds", 
                              chains = 2, cores = 2, iter = 2000,
                              thin = 4, ...){
  library(rstan) # prevent annoying NAMESPACE issue described in
  # https://github.com/stan-dev/rstan/issues/353
  
  set.seed(seed)
  settings = yaml::yaml.load_file("settings.yaml")
  
  # If we don't have any observed richness, we can't use the row, so omit
  # anything with is.na(richness)
  df = get_pop_ts_env_data(settings$start_yr,
                           settings$end_yr,
                           settings$last_train_year,
                           settings$min_year_percentage) %>% 
    collapse_to_richness() %>% 
    filter(!is.na(richness)) %>%
    mutate(in_train = year <= settings$last_train_year) %>% 
    add_observer_index() %>% 
    add_site_index()
  
  # Stan only gets the training data, plus some counts (N, N_observer, etc.)
  stan_data = df %>% 
    filter(in_train) %>% 
    select(observer_index, site_index, richness) %>% 
    as.list() %>% 
    c(N = sum(df$in_train), 
      N_observer = n_distinct(.$observer_index),
      N_site = n_distinct(.$site_index)) %>% 
    c(N_test_observer = n_distinct(df$observer_index) - n_distinct(.$observer_index))
  
  # Compile the model
  model = rstan::stan_model(stan_file)
  
  # Fit the model; drop character vectors that Stan doesn't need/like
  samples = rstan::sampling(model, 
                            data = stan_data, 
                            chains = chains, cores = cores, iter = iter,
                            thin = thin,
                            ...)
  
  # Save tidy tables of observer effects and site effects
  observer_table = tidy_stan(samples, 
                             c("observer_effect", "test_observer_effect"), 
                             key = "observer_index", 
                             value = "observer_effect")
  site_table = tidy_stan(samples, 
                         "site_effect", 
                         key = "site_index", 
                         value = "site_effect")
  
  
  # Store the model output
  obs_model = c(
    rstan::extract(samples),
    data = list(
      df %>% 
        left_join(observer_table, "observer_index") %>% 
        left_join(site_table, c("site_index", "iteration"))
    )
  )
  # Print the R-hat MCMC diagnostic
  max_rhat = max(rstan::summary(samples)[[1]][,"Rhat"])
  message("R-hat statistics all at or below ", round(max_rhat, 2))
  message("R-hat values should be close to 1; values larger than about 1.05\nindicate convergence problems")
  
  saveRDS(obs_model, file = output_file)
  saveRDS(samples, file = paste0("rstan_", output_file))
  
  NULL
}

#####
# Helper functions
#####

# Make a tidy table out of one or more parameter matrices
tidy_stan = function(samples, names, key = key, value = value){
  # Extract the matrices, convert them to data_frames, and bind them into one
  # data_frame.  The column names correspond to index values (1:N)
  df = rstan::extract(samples, names) %>% 
    purrr::map(dplyr::as_data_frame) %>% 
    dplyr::bind_cols() %>% 
    purrr::set_names(1:ncol(.))
  
  # Initially, iteration number is denoted by row; make this explicit and then
  # make the data tidy/long instead of wide.
  out = df %>% 
    cbind(iteration = 1:nrow(.)) %>% 
    gather(key = "key", value = "value", -iteration)
  
  # the "key" column (i.e. the index) needs to be coerced from character because
  # it was stored as a column name before `gather`ing.
  # Nonstandard evaluation makes it hard to pass key/value names in to `gather`,
  # so I add them here.
  out %>% 
    mutate(key = as.integer(key)) %>% 
    purrr::set_names(c("iteration", key, value))
}


add_site_index = function(x){
  # index the sites by their site_id from 1 to n_sites.
  x %>% 
    distinct(site_id) %>% 
    mutate(site_index = 1:nrow(.)) %>% 
    right_join(x, "site_id")
}


add_observer_index = function(x) {
  # The Stan model assumes that observers are indexed by sequential whole numbers
  # and that the observers that participated in the training set come first in
  # the sequence.
  x %>% 
    group_by(observer_id) %>% 
    summarize(training_observer = any(in_train)) %>% 
    arrange(dplyr::desc(training_observer)) %>% 
    mutate(observer_index = 1:nrow(.)) %>% 
    right_join(x, "observer_id")
}
