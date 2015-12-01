library(gbm)
library(parallel)
set.seed(1)

mc.cores = 8

train_x = as.data.frame(readRDS("data/train_x.rds"))
validation_x = as.data.frame(readRDS("data/validation_x.rds"))
train_y = readRDS("data/train_y.rds")

predictions = mclapply(
  1:ncol(train_y),
  function(i){
    g = gbm(
      formula = train_y[ , i] ~ ., 
      distribution = "poisson", 
      data = train_x,
      n.trees = 1E4,
      interaction.depth = 3,
      cv.folds = 5,
      shrinkage = .001,
      verbose = TRUE,
      n.cores = 1
    )
    
    perf = gbm.perf(g, method = "cv")
    predict(g, validation_x, n.trees = perf, type = "response")
  },
  mc.cores = mc.cores
)

write.csv(predictions, file = "gbm_predictions.csv")
