print(mistnet::adam.updater()$annealing_rate) # will fail if mistnet version is old

library(doParallel)
cl=makeCluster(16)
registerDoParallel(cl)

foreach(i=1:500) %dopar% {
  source("mistnet.R")
  fit_mistnet(iter = i, use_obs_model = TRUE)
}

foreach(i=1:500) %dopar% {
  source("mistnet.R")
  fit_mistnet(iter = i, use_obs_model = FALSE)
}
