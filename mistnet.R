library(Heraclitus)
library(tidyverse)
library(mistnet)
library(progress)
devtools::load_all() 


timeframe = "train_32"
prepend_timeframe = function(x) {
  paste0("results", "/", timeframe, "/", x)
}


# Ingest the `settings` and put the timeframe of interest at the top level
# of the list.
settings = yaml::yaml.load_file("settings.yaml")
settings = c(settings, settings$timeframes[[timeframe]])
settings$timeframes = NULL
is_future = settings$timeframe == "future"

bbs = get_bbs_data() %>% 
  mutate(species_id = paste0("sp", species_id), presence = abundance > 0) %>% 
  select(-lat, -long, -abundance) %>% 
  spread(key = species_id, value = presence, fill = 0)

# Mistnet settings based on the second-best model in Appendix D of 
# the mistnet paper

# Other mistnet settings
n_hours = 12
iteration_size = 10
n_prediction_samples = 500

# If CV==TRUE, we're cross-validating within the training set.
# a_0 and annealing_rate are hyperparameters of the `adam` optimizer
fit_mistnet = function(iter,
                       use_obs_model, 
                       mistnet_arglist,
                       updater_arglist,
                       CV = FALSE){
  
  # Collect data for this iteration, then discard the rest to save
  # memory
  obs_model = readRDS(prepend_timeframe("observer_model.rds"))
  data = obs_model$data %>% 
    filter(iteration == iter) %>% 
    left_join(bbs, c("site_id", "year"))
  
  if (is_future) {
    future = get_env_data(timeframe = "future") %>% 
      filter(site_id %in% !!data$site_id,
             year > !!settings$last_train_year) %>% 
      select(which(colnames(.) %in% !!colnames(data))) %>% 
      mutate(observer_effect = rnorm(nrow(.), mean = 0, 
                                     sd = obs_model$observer_sigma[[iter]]),
             in_train = FALSE)
    data = bind_rows(data, future)
    
    rm(future)
  }
  rm(obs_model)
  gc(TRUE)
  
  
  
  # Cross-validation
  if (CV) {
    # We only get to see the data up through the last_train_year
    data = filter(data, year <= settings$last_train_year)
    
    # All years until the last_train_year are in the training set for CV
    data$in_train = data$year != settings$last_train_year
  }
  
  vars = c(settings$vars, if (use_obs_model) {"observer_effect"})
  
  x = data %>% 
    select(one_of(vars)) %>% 
    as.matrix()
  
  # Z-scale each column based on mean & sd of training set
  for (i in 1:ncol(x)) {
    mean_i = mean(x[data$in_train, i])
    sd_i = sd(x[data$in_train, i])
    x[,i] = (x[,i] - mean_i) / sd_i
  }
   
  # response variables start wtih "sp" followed by numbers & occur at least 
  # once.
  y = data %>% 
    select(matches("^sp[0-9]+$")) %>% 
    select(which(colSums(. , na.rm = TRUE) > 0)) %>% 
    as.matrix()
  
  # Drop columns we won't need below
  data = data %>% 
    select(site_id, year, iteration, richness, in_train)
  gc(TRUE)
  
  # Model definition --------------------------------------------------------
  
  net = mistnet(
    x = x[data$in_train, ],
    y = y[data$in_train, ],
    layer.definitions = list(
      defineLayer(
        nonlinearity = elu.nonlinearity(),
        size = mistnet_arglist$N1,
        prior = gaussian.prior(mean = 0, sd = .5)
      ),
      defineLayer(
        nonlinearity = elu.nonlinearity(),
        size = mistnet_arglist$N2,
        prior = gaussian.prior(mean = 0, sd = .5)
      ),
      defineLayer(
        nonlinearity = sigmoid.nonlinearity(),
        size = ncol(y),
        prior = gaussian.prior(mean = 0, sd = .5)
      )
    ),
    loss = bernoulliRegLoss(a = 1 + 1E-6, b = 1 + 1E-6),
    updater = purrr::invoke(adam.updater$new, updater_arglist),  
    sampler = gaussian.sampler(ncol = mistnet_arglist$latent_dim, sd = 1),
    n.importance.samples = mistnet_arglist$n.importance.samples,
    n.minibatch = mistnet_arglist$n.minibatch,
    training.iterations = 0,
    initialize.biases = TRUE,
    initialize.weights = TRUE
  )
  
  # Fit the model ---------------------------------------------------------
  print("entering training loop")
  start_time = Sys.time()
  while (difftime(Sys.time(), start_time, units = "hours") < n_hours) {
    net$fit(iteration_size)
    # Update prior variance
    for (layer in net$layers) {
      layer$prior$update(
        layer$weights, 
        update.mean = FALSE, 
        update.sd = TRUE,
        min.sd = .01
      )
      layer$prior$sd = layer$prior$sd * mistnet_arglist$sd_mult
    }
    # Update mean for final layer
    net$layers[[3]]$prior$update(
      layer$weights, 
      update.mean = TRUE, 
      update.sd = FALSE,
      min.sd = .01
    )
    if (any(is.na(net$layers[[3]]$outputs))) {
      stop("NA in network output")
    }
  } # End while
  
  print("completed training loop")
  
  gc(TRUE)
  
  # Streaming mean & variance in expected values
  moments = moment_stream$new() 
  residual_variance = 0 # Variance in richness, given occurrence probabilities
  ll = matrix(NA, sum(!data$in_train, na.rm = TRUE), n_prediction_samples)
  newdata = as.matrix(x[!data$in_train, ])

  for (i in 1:n_prediction_samples) {
    
    p = predict(net, 
                newdata = newdata, 
                n.importance.samples = 1)
    
    if (!is_future) {
      ll[,i] = rowSums(dbinom(x = y[!data$in_train, ], 
                              size = 1, 
                              prob = p, 
                              log = TRUE))
    }
    
    moments$update(list(rowSums(p)))
    residual_variance = residual_variance + 
      rowSums(p * (1 - p)) / n_prediction_samples
  }
  print("predictions generated")
  
  # total variance is variance in mean estimates plus mean of the variance 
  # estimates
  out = data %>% 
    filter(!in_train) %>% 
    cbind(mean = moments$m,
          sd = sqrt(moments$v_hat + residual_variance)) %>% 
    select(-in_train) %>% 
    mutate(model = "mistnet", use_obs_model = use_obs_model,
           log_lik = apply(ll, 1, mistnet:::logMeanExp))
  
  # Save the predictions
  dir.create("mistnet_output", showWarnings = FALSE)
  saveRDS(out, file = prepend_timeframe(paste0("mistnet_output/", "iteration_", 
                                               iter, "_", 
                                               ifelse(CV, "CV", use_obs_model), 
                                               ".rds")))
}
