#!/usr/bin/Rscript
library(tidyverse)
source("mistnet.R")
args <- commandArgs(TRUE)
start = args[[1]]
end = args[[2]]
timeframe = args[[3]]

# Find the cross-validation iteration with the lowest error
row = dir(paste0("results/", timeframe, "/mistnet_output/"), 
          pattern = "CV", 
          full.names = TRUE) %>% 
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

for (i in start:end) {
  for (use_obs_model in c(TRUE, FALSE)) {
    gc(TRUE)
    fit_mistnet(iter = i, use_obs_model = use_obs_model, CV = FALSE,
                mistnet_arglist = hyper, updater_arglist = hyper)
  }
}
