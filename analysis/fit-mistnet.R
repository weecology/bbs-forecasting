source("mistnet.R")
library(parallel)
mclapply(1:500,
         function(i){
           fit_mistnet(iter = i, use_obs_model = TRUE)
         },
         mc.cores = 16)
