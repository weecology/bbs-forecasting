# Install pacman if it isn't already installed

if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")


# Install analysis packages using pacman

pacman::p_load(dplyr, forecast, ggplot2, Hmisc, tidyr, mgcv, sp, raster,
               maptools, doParallel, stringr, RCurl, roxygen2,
	       broom, devtools, doParallel, dplyr, forecast, ggplot2,
	       gimms, Hmisc, maptools, mgcv, prism, raster, stringr, sp,
	       tidyr, rgdal, rgeos, DBI, RSQLite, lme4, caret, mapproj,
               viridis, git2r, rstan, readr, purr, gbm, randomForest)
pacman::p_load_gh('seantuck12/MODISTools', 'ropensci/ecoretriever')


# Install IRkernal packages and spec

pacman::p_load(repr, pbdZMQ)
pacman::p_load_gh('IRkernel/IRdisplay', 'IRkernel/IRkernel')
IRkernel::installspec()
