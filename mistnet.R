library(Heraclitus)
library(tidyverse)
library(mistnet)
library(progress)
devtools::load_all()
if (!file.exists("observer_model.rds")) {
  fit_observer_model()
}

bbs = get_bbs_data() %>% 
  mutate(species_id = paste0("sp", species_id), presence = abundance > 0) %>% 
  select(-lat, -long, -abundance) %>% 
  spread(key = species_id, value = presence, fill = 0)

# Mistnet settings based on the second-best model in Appendix D of 
# the mistnet paper
n.minibatch = 79L
latent_dim = 10L
n.importance.samples = 28L
N1 = 37L
N2 = 15L

# Other mistnet settings
n_opt_iterations = 1E4
iteration_size = 10
n_prediction_samples = 500

# If CV==TRUE, we're cross-validating within the training set.
# a_0 and annealing_rate are hyperparameters of the `adam` optimizer
fit_mistnet = function(iter,
                       use_obs_model, 
                       updater_arglist,
                       CV = FALSE){
  
  settings = yaml::yaml.load_file("settings.yaml")
  
  # Collect data for this iteration, then discard the rest to save
  # memory
  obs_model = readRDS("observer_model.rds")
  data = filter(obs_model$data, iteration == iter) %>% 
    left_join(bbs, c("site_id", "year"))
  rm(obs_model)
  gc()
  
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
  
  y = data %>% 
    select(matches("^sp[0-9]+$")) %>% 
    select(which(colSums(.) > 0)) %>% 
    as.matrix()
  
  # Drop columns we won't need below
  data = data %>% 
    select(site_id, year, iteration, richness, in_train)
  gc()
  
  # Model definition --------------------------------------------------------
  
  net = mistnet(
    x = x[data$in_train, ],
    y = y[data$in_train, ],
    layer.definitions = list(
      defineLayer(
        nonlinearity = rectify.nonlinearity(),
        size = N1,
        prior = gaussian.prior(mean = 0, sd = .5)
      ),
      defineLayer(
        nonlinearity = linear.nonlinearity(),
        size = N2,
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
    sampler = gaussian.sampler(ncol = latent_dim, sd = 1),
    n.importance.samples = n.importance.samples,
    n.minibatch = n.minibatch,
    training.iterations = 0,
    initialize.biases = TRUE,
    initialize.weights = TRUE
  )
  # First layer biases equal 1; prevents "dead" rectifiers
  net$layers[[1]]$biases[] = 1
  
  # Shrink each layer's initial weights compared to baseline
  net$layers[[1]]$weights = net$layers[[1]]$weights / 2
  net$layers[[2]]$weights = net$layers[[2]]$weights / 2
  net$layers[[3]]$weights = net$layers[[3]]$weights / 2
  
  # Fit the model ---------------------------------------------------------
  pb = progress_bar$new(format = "[:bar] :percent eta: :eta",
                        total = n_opt_iterations / iteration_size)
  pb$tick(0)
  while (net$completed.iterations < n_opt_iterations) {
    net$fit(iteration_size)
    pb$tick()
    # Update prior variance
    for (layer in net$layers) {
      layer$prior$update(
        layer$weights, 
        update.mean = FALSE, 
        update.sd = TRUE,
        min.sd = .01
      )
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
  
  # Streaming mean & variance in expected values
  moments = moment_stream$new() 
  residual_variance = 0 # Variance in richness, given occurrence probabilities
  for (i in 1:n_prediction_samples) {
    p = predict(net, 
                newdata = x[!data$in_train, ], 
                n.importance.samples = 1)
    
    moments$update(list(rowSums(p)))
    residual_variance = residual_variance + 
      rowSums(p * (1 - p)) / n_prediction_samples
  }
  
  # total variance is variance in mean estimates plus mean of the variance 
  # estimates
  out = data %>% 
    filter(!in_train) %>% 
    cbind(mean = moments$m,
          sd = sqrt(moments$v_hat + residual_variance)) %>% 
    select(-in_train) %>% 
    mutate(model = "mistnet", use_obs_model = use_obs_model)
  
  dir.create("mistnet_output", showWarnings = FALSE)
  saveRDS(out, file = paste0("mistnet_output/", "iteration_", iter, 
                             "_", ifelse(CV, "CV", use_obs_model), ".rds"))
}

# # Plotting species' responses to environmental variables
# N = 250
# xy = as.matrix(expand.grid(seq(-4, 4, length = N), seq(-4, 4, length = N)))
# xx = cbind(as.matrix(x[rep(1, N^2), ]), matrix(0, N^2, latent_dim))
# xx[, c(3, 11)] = xy
# 
# relu = function(x){
#   structure(pmax(x, 0), dim = dim(x))
# }
# 
# to_plot = xx %*% 
#   net$layers[[1]]$weights %>% 
#   relu() %*% 
#   net$layers[[2]]$weights %*% net$layers[[3]]$weights %>% 
#   plogis()
# 
# cbind(xy, fill = to_plot[,28]) %>% 
#   as.data.frame() %>% 
#   ggplot(aes(x = Var1, y = Var2, fill = fill)) + 
#   geom_raster() +
#   viridis::scale_fill_viridis(limits = c(0, 1)) + 
#   cowplot::theme_cowplot()

# Summarizing mean squared error
# out %>%
#   group_by(year, site_id, richness) %>%
#   summarize(mean = mean(mean)) %>%
#   summarize(mse = mean((richness - mean)^2)) %>%
#   group_by(year) %>%
#   summarize(mse = mean(mse))

