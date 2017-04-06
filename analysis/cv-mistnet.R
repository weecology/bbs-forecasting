#!/usr/bin/Rscript
source("mistnet.R")
args <- commandArgs(TRUE)
start = args[[1]]
end = args[[2]]
N = args[[3]]

discrete_log_runif = function(N, a, b){
  # Sample from a discrete log-uniform distribution between a and b
  # Note that as.integer rounds down, so the upper limit needs to be
  # incremented by 1.
  as.integer(exp(runif(N, log(a), log(b + 1))))
}

# check for mistnet version: this will fail if mistnet version is too old.
# (Otherwise, bogus arguments will get absorbed into `...` without warning)
invisible(mistnet::adam.updater()$annealing_rate)

# Store the current random seed so we don't permanently overwrite it
seed = .Random.seed

# Deterministically set hyperparameter list by setting the seed
set.seed(1)
mistnet_arglist = list(n.minibatch = discrete_log_runif(N, 10, 100),
                       latent_dim =  discrete_log_runif(N, 5, 20),
                       n.importance.samples = discrete_log_runif(N, 10, 50),
                       N1 = discrete_log_runif(N, 40, 100),
                       N2 = discrete_log_runif(N, 20, 50),
                       sd_mult = rgamma(N, 10, 5))

updater_arglist = list(a_0 = rlnorm(N, log(.005), 1),
                    annealing_rate = rlnorm(N, -10, 2),
                    b1 = rbeta(N, 12, 2),
                    b2 = rbeta(N, 100, 1),
                    e = rlnorm(N, -18, 2))

# Re-set the seed to its previous value
set.seed(seed)

for (i in start:end) {
  gc(TRUE)
  # When cross-validating, don't use the observation model because that would
  # add noise across iterations with different MCMC estimates of the observer
  # effects.
  fit_mistnet(iter = i, 
              use_obs_model = FALSE, 
              mistnet_arglist = purrr::map(mistnet_arglist, i),
              updater_arglist = purrr::map(updater_arglist, i), 
              CV = TRUE)
}
