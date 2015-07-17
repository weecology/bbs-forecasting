# Get NDVI data from MODIS

library(MODISTools)
library(dplyr)

route_data <- read.csv('data/BBS_routes.csv')
route_data$site_id <-  route_data$statenum * 1000 + route_data$route

contig_modern_routes <- read.csv('data/contig_modern_routes.csv')

colnames(contig_modern_routes) <- c('site_id')
route_data <- semi_join(route_data, contig_modern_routes)


coord_data <- route_data[c('lati', 'loni')]
coord_data <- unique(coord_data) # the Routes table has duplicate coordinate values for some reason
colnames(coord_data) <- c('lat', 'long')
coord_data$start.date <- rep(2000, nrow(coord_data))
coord_data$end.date <- rep(2014, nrow(coord_data))

unaquired_coord_data <- UpdateSubsets(LoadDat = coord_data, StartDate = TRUE,
                                     Dir = "./data/modisdata/")

MODISSubsets(LoadDat = unaquired_coord_data, Products = "MOD13Q1",
             Bands = c("250m_16_days_NDVI"), Size = c(1,1),
             SaveDir = "./data/modisdata", StartDate = TRUE)
