library(tidyverse)
library(mistnet)
library(progress)
devtools::load_all()
settings = yaml::yaml.load_file("settings.yaml")
n_prediction_samples = 20
use_obs_model = TRUE

if (!file.exists("observer_model.rds")) {
  fit_observer_model()
}
obs_model = readRDS("observer_model.rds")

bbs = get_bbs_data() %>% 
  mutate(species_id = paste0("sp", species_id), presence = abundance > 0) %>% 
  select(-lat, -long, -abundance) %>% 
  spread(key = species_id, value = presence, fill = 0)


i = 1
data = filter(obs_model$data, iteration == i) %>% 
  left_join(bbs, c("site_id", "year"))

# Free up some memory
rm(bbs, obs_model);gc()

vars = yaml::yaml.load_file("settings.yaml")$vars %>% 
  c(if (use_obs_model) {"observer_effect"})

x = data %>% 
  select(one_of(vars)) %>% 
  as.matrix()
# Scaling incorrectly because it's quicker for now
x = scale(x)

y = data %>% 
  select(matches("^sp[0-9]+$")) %>% 
  select(which(colSums(.) > 0)) %>% 
  as.matrix()



# Model definition --------------------------------------------------------
N1 = 45L
N2 = 15L
rate = 0.1
latent_dim = 6L
n.importance.samples = 30L
n.minibatch = 25L

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
  updater = adam.updater(a = 0.01, b2 = 0.99),
  sampler = gaussian.sampler(ncol = latent_dim, sd = 1),
  n.importance.samples = n.importance.samples,
  n.minibatch = n.minibatch,
  training.iterations = 0,
  initialize.biases = TRUE,
  initialize.weights = TRUE
)
# First layer biases equal 1; prevents "dead" rectifiers
net$layers[[1]]$biases[] = 1

# Shrink each layer's initail weights
net$layers[[1]]$weights = net$layers[[1]]$weights / 2
net$layers[[2]]$weights = net$layers[[2]]$weights / 2
net$layers[[3]]$weights = net$layers[[3]]$weights / 2

# save weights
initial_weights = net$layers %>% map(~.x[["weights"]])

# Fit the model -----------------------------------------------------------
n_opt_iterations = 1E3
iteration_size = 10
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
  ll = sapply(1:n.importance.samples,
              function(i){
                net$loss$loss(y = net$y[net$row.selector$minibatch.ids, ], 
                              yhat = net$layers[[3]]$outputs[,,i]) %>% 
                  rowSums()
              })
  print(mean(ll * net$importance.weights))
} # End while


hist(apply(net$layers[[3]]$outputs, 2, rowMeans), breaks = "fd", border = 2, yaxs = "i")
hist(qlogis(apply(net$layers[[3]]$outputs, 2, rowMeans)), breaks = "fd", border = 2, yaxs = "i")
hist(net$layers[[3]]$coef.updater$delta, breaks = "fd", border = 2)

p = predict(net, newdata = x[!data$in_train, ], n.importance.samples = n_prediction_samples)

make_mistnet_tidy = function(p, value_name){
  apply(p, 1, colSums) %>% 
    t() %>% 
    cbind(data[!data$in_train, c("site_id", "year", "iteration")]) %>% 
    gather_(key_col = "sample", value_col = value_name, 
            gather_cols = as.character(1:n_prediction_samples))
}

out = make_mistnet_tidy(p, "mean") %>% 
  left_join(make_mistnet_tidy(p * (1 - p), "variance"),
            c("site_id", "year", "iteration", "sample")) %>% 
  left_join(select(data, site_id, year, richness), c("site_id", "year"))

out %>% 
  group_by(year, site_id, richness) %>% 
  summarize(mean = mean(mean)) %>% 
  summarize(mse = mean((richness - mean)^2)) %>%  
  group_by(year) %>% 
  summarize(mse = mean(mse))
