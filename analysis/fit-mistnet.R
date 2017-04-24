#!/usr/bin/Rscript
library(tidyverse)
source("mistnet.R")
args <- commandArgs(TRUE)
start = args[[1]]
end = args[[2]]

# Find the cross-validation iteration with the lowest error
row = dir("mistnet_output/", pattern = "CV", full.names = TRUE) %>% 
  map(readRDS) %>% 
  bind_rows() %>% 
  group_by(iteration) %>% 
  summarize(mse = Metrics::mse(richness, mean)) %>% 
  .[["mse"]] %>% 
  which.min()
  

# Grab the hyperparameters that produced the lowest error
hyper = read.csv("mistnet_hyper.csv") %>% 
  filter(X == row) %>% 
  as.list()

fit_mistnet(iter = i, 
            use_obs_model = FALSE, 
            mistnet_arglist = purrr::map(mistnet_arglist, i),
            updater_arglist = purrr::map(updater_arglist, i), 
            CV = TRUE)

for (i in start:end) {
  for (use_obs_model in c(TRUE, FALSE)) {
    gc(TRUE)
    fit_mistnet(iter = i, use_obs_model = use_obs_model, CV = FALSE,
                mistnet_arglist = hyper, updater_arglist = hyper)
  }
}
