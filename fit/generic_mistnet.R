full_bbs = na.omit(readRDS("data/munged.Rds"))
library(mistnet)
library(fastICA)
library(dplyr)

focal_year = 1990

set.seed(1)

# Train on the training data through focal_year
train_bbs = full_bbs[full_bbs$in_train & full_bbs$year <= focal_year, ]

make_rescaler = function(x){
  scaled_x = scale(x)
  
  # Half the variance in X should come from temp and half from precip,
  # but 3/4 of the columns are about temperature.
  temp_columns = grep("^t", colnames(scaled_x))
  scaled_x[ , temp_columns] = scaled_x[ , temp_columns] / sqrt(3)
  
  ica = fastICA(scaled_x, 10)
  
  function(new_x){
    # Z-transform the new variables using the old means and sds
    scaled_new_x = sapply(
      1:ncol(new_x), 
      function(i){
        (new_x[[i]] - attr(scaled_x,"scaled:center")[i]) / attr(scaled_x, "scaled:scale")[i]
      }
    )

    # Shrink the temperature-related variables as above
    scaled_new_x[ , temp_columns] = scaled_new_x[ , temp_columns] / sqrt(3)

    # Project the new x variables according to the ICA built on the training data
    scaled_new_x %*% ica$K %*% ica$W
  }
}


# Validate on the test data through focal_year
validation_bbs = full_bbs[full_bbs$in_test & full_bbs$year <= focal_year, ]

# Once we have hyperparameters, fit to everything through focal year
fitting_bbs = full_bbs[full_bbs$year <= focal_year, ]

# Finally, test on all the post-focal_year data
test_bbs = full_bbs[full_bbs$year > focal_year, ]


species_ids = tbl(src_sqlite("data/bbsforecasting.sqlite"), "bbs_species") %>% 
  collect() %>% 
  extract2("species_id")

climate_colnames = grep("ppt|tmin|tmean|tmax", colnames(train_bbs), value = TRUE)
species_colnames = colnames(full_bbs)[colnames(full_bbs) %in% species_ids]


train_rescaler = make_rescaler(train_bbs[ , climate_colnames])

train_x = train_rescaler(train_bbs[ , climate_colnames])
validation_x = train_rescaler(validation_bbs[ , climate_colnames]) 


train_y = as.matrix(train_bbs[ , species_colnames])
train_y = train_y[ , colSums(train_y) > 0]

validation_y = as.matrix(validation_bbs[ , colnames(train_y)])

# Hyperparameters ---------------------------------------------------------

# Roughly based on what worked well in the Harris 2015 paper
size1 = 35L        # Number of units in the first hidden layer
size2 = 12L        # Number of units in the second hidden layer
n_samples = 30L    # Number of importance samples
rate = 0.1         # Learning rate 
n_minibatch = 75L  # Number of routes in a minibatch
n_latent = 10L     # Number of latent variables


net = mistnet(
  x = train_x,
  y = train_y,
  layer.definitions = list(
    defineLayer(
      nonlinearity = rectify.nonlinearity(),
      size = size1,
      prior = gaussian.prior(mean = 0, sd = .1)
    ),
    defineLayer(
      nonlinearity = linear.nonlinearity(),
      size = size2,
      prior = gaussian.prior(mean = 0, sd = .1)
    ),
    defineLayer(
      nonlinearity = exp.nonlinearity(),
      size = ncol(train_y),
      prior = gaussian.prior(mean = 0, sd = .1)
    )
  ),
  loss = poissonLoss(),
  updater = adagrad.updater(learning.rate = rate),
  sampler = gaussian.sampler(ncol = n_latent, sd = 1),
  n.importance.samples = n_samples,
  n.minibatch = n_minibatch,
  training.iterations = 0,
  initialize.biases = TRUE,
  initialize.weights = TRUE
)

# Set first layer biases to 1 to prevent "dead" RELUs
net$layers[[1]]$biases[] = 1

# Train the model
for(i in 1:100){
  net$fit(10) # Climb the likelihoood surface
  cat(".")
  
  # Update prior variance in each layer
  for(layer in net$layers){
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
}


plot(
  t(apply(net$layers[[3]]$outputs, 1, rowMeans)),
  net$y[net$row.selector$minibatch.ids, ],
  log = "x"
)
x_seq = 10^seq(-5, 5, length = 128)
lines(x_seq, x_seq, col = 2)

