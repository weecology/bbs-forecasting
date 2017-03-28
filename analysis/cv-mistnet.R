#!/usr/bin/Rscript
source("mistnet.R")
args <- commandArgs(TRUE)
start = args[[1]]
end = args[[2]]
N = args[[3]]

# check for mistnet version: this will fail if mistnet version is too old.
# (Otherwise, bogus arguments will get absorbed into `...` without warning)
invisible(mistnet::adam.updater()$annealing_rate)

# Store the current random seed so we don't permanently overwrite it
seed = .Random.seed

# Deterministically set hyperparameter list by setting the seed
set.seed(1)
full_arglist = list(a_0 = rlnorm(N, -5, 1),
                    annealing_rate = rlnorm(N, -7, 1.5),
                    b1 = rbeta(N, 5, 1),
                    b2 = rbeta(N, 100, 1),
                    e = rlnorm(N, -10, 4))
# Re-set the seed to its previous value
set.seed(seed)

for (i in start:end) {
  # When cross-validating, don't use the observation model because that would
  # add noise across iterations with different MCMC estimates of the observer
  # effects.
  fit_mistnet(iter = i, 
              use_obs_model = FALSE, 
              updater_arglist = purrr::map(full_arglist, i), CV = TRUE)
}
