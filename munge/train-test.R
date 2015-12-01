full_bbs = na.omit(readRDS("data/munged.Rds"))
library(fastICA)
library(magrittr)
library(dplyr)
library(beepr)

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
train_y = train_y[ , colSums(train_y) > 0] # Drop empty columns

validation_y = as.matrix(validation_bbs[ , colnames(train_y)])

saveRDS(train_x, file = "data/train_x.rds")
saveRDS(validation_x, file = "data/validation_x.rds")
saveRDS(train_y, file = "data/train_y.rds")
saveRDS(validation_y, file = "data/validation_y.rds")

