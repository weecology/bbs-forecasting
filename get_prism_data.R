# Download and install PRISM data into PostgreSQL and values for BBS routes

library(prism)
library(raster)
library(dplyr)
library(tidyr)
options(prism.path = "./data/prismdata_retriever")

months <- c(1:12)
clim_vars <- c("ppt", "tmin", "tmean", "tmax")
for (month in months){
  for (clim_var in clim_vars){
    get_prism_monthlys(type=clim_var, year = 1981:2014, month = month, keepZip=F)
  }
}

datadirs = dir(datapath)
for (datadir in datadirs) {
  bil_file = paste(datadir, '.bil', sep = "")
  bil_file_path = file.path(datapath, datadir, bil_file)
  sql_file = paste(datadir, ".sql", sep = "")
  sql_file_path = file.path(datapath, datadir, sql_file)
  system(paste("raster2pgsql -s 4326", bil_file_path, datadir, ">", sql_file_path))
  system(paste("psql -d bbsforecasting -f", sql_file_path))
}

bbs_data <- read.csv("data/bbs_data.csv")
locations <- unique(select(bbs_data, site_id, long, lat))
coordinates(locations) <- c("long", "lat")
#prism_files <- grep("_[0-9]{4}[0]",ls_prism_data()[,1],value=T)
prism_stacked <- prism_stack(ls_prism_data())
extracted <- raster::extract(prism_stacked, locations)
prism_bbs_data <- data.frame(site_id = locations$site_id, coordinates(locations), extracted)
prism_bbs_data <- prism_bbs_data %>%
                    gather(date, value, 4:ncol(prism_bbs_data)) %>%
                    tidyr::extract(date, c("clim_var", "year", "month"),
                                   "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_") %>%
                    spread(clim_var, value)
prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
prism_bbs_data$month <- as.numeric(prism_bbs_data$month)
database <- src_sqlite("./data/bbsforecasting.sqlite", create = TRUE)
mydata <- copy_to(database, prism_bbs_data, temporary = FALSE,
                  indexes = list(c("site_id", "year", "month")))
