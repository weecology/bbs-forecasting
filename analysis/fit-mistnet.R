#!/usr/bin/Rscript
source("mistnet.R")
args <- commandArgs(TRUE)
start = args[[1]]
end = args[[2]]

# check for mistnet version: this will fail if mistnet version is too old
# (Otherwise, bogus arguments will get absorbed into `...` without warning)
invisible(mistnet::adam.updater()$annealing_rate)

for (i in start:end) {
  gc(TRUE)
  fit_mistnet(iter = i, use_obs_model = TRUE)
  gc(TRUE)
  fit_mistnet(iter = i, use_obs_model = FALSE)
}
