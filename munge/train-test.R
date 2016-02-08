library(fastICA)
library(magrittr)
library(dplyr)

full_bbs = na.omit(readRDS("data/munged.rds"))
focal_year = readRDS("data/focal_year.rds")

set.seed(1)

# Split rows --------------------------------------------------------------

# Train on the training data through focal_year
train_bbs = full_bbs[full_bbs$in_train & full_bbs$year <= focal_year, ]

# Validate on the test data through focal_year
validation_bbs = full_bbs[full_bbs$in_test & full_bbs$year <= focal_year, ]

# Once we have hyperparameters, fit to everything through focal year
fitting_bbs = full_bbs[full_bbs$year <= focal_year, ]

# Finally, test on all the post-focal_year data
test_bbs = full_bbs[full_bbs$year > focal_year, ]


# Identify columns --------------------------------------------------------

species_ids = tbl(src_sqlite("data/bbsforecasting.sqlite"), "bbs_species") %>% 
  collect() %>% 
  extract2("species_id")

climate_colnames = grep("ppt|tmin|tmean|tmax", colnames(train_bbs), value = TRUE)
species_colnames = colnames(full_bbs)[colnames(full_bbs) %in% species_ids]

is_temp = grepl("^t", climate_colnames)


# Environmental munging functions -----------------------------------------

# Like `scale`, but dividing by (sd * sd_factor) instead of just by sd
env_scale = function(x, sd_factor = 1) {
  out = scale(x)
  
  for (i in 1:ncol(x)) {
    out[ , i] = out[ , i] / sd_factor[i]
  }
  attr(out, "scaled:scale") = attr(out, "scaled:scale") * sd_factor
  
  out
}

# Adust the out_of_sample values by subtracting off the mean and dividing by the
# standard deviation used to scale the pre_scaled values.
rescale = function(pre_scaled, out_of_sample) {
  means = attr(pre_scaled, "scaled:center")
  sds = attr(pre_scaled, "scaled:scale")
  
  # Convert for safety because data.frames, matrices, and tbl_dfs can respond
  # differently to `[`
  out_of_sample = as.matrix(out_of_sample)
  
  out = sapply(
    1:ncol(pre_scaled),
    function(i) {
      (out_of_sample[ , i] - means[i]) / sds[i]
    }
  )
  colnames(out) = colnames(out_of_sample)
  
  out
}

# Project (rotate) the scaled x variables according to the ICA model
project_ica = function(scaled_x, ica_model) {
  scaled_x %*% ica_model$K %*% ica_model$W
}


# Rescale and rotate ------------------------------------------------------

# Will eventually repeat this process for fitting and test data, but it's 
# probably better not to calculate them until we need them so we don't 
# accidentally peek.

# Dimensionality to use for ICA
n_components = 10

# Three times as many columns are temperature-related as precip related. 
# Dividing the precip variance by three (i.e. reducing the sd by sqrt(3)) makes
# temp have equal variance to precip
sd_factor = ifelse(is_temp, sqrt(3), 1)

# Scaled (but not rotated)
scaled_train_x = env_scale(train_bbs[ , climate_colnames], sd_factor = sd_factor)
scaled_validation_x = rescale(
  pre_scaled = scaled_train_x,
  out_of_sample = validation_bbs[ , climate_colnames]
)

# Optimize the rotation on the scaled training data
train_ica = fastICA(scaled_train_x, n_components)

# Apply the same rotation to the training and validation data
train_x = project_ica(scaled_train_x, train_ica)
validation_x = project_ica(scaled_validation_x, train_ica) 


# Extract response variables ----------------------------------------------

train_y = as.matrix(train_bbs[ , species_colnames])
train_y = train_y[ , colSums(train_y) > 0] # Drop empty columns

validation_y = as.matrix(validation_bbs[ , colnames(train_y)])


# Save the outputs --------------------------------------------------------

saveRDS(train_x, file = "data/train_x.rds")
saveRDS(validation_x, file = "data/validation_x.rds")
saveRDS(train_y, file = "data/train_y.rds")
saveRDS(validation_y, file = "data/validation_y.rds")

