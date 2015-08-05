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
coord_date_data <- coord_data
coord_date_data$start.date <- rep(as.Date("2000-06-01"), nrow(coord_data))
coord_date_data$end.date <- rep(as.Date("2000-06-30"), nrow(coord_data))

years = 2001:2015
for (year in years){
  coord_data_year <- coord_data
  coord_data_year$start.date <- rep(as.Date(paste(year, "-06-01", sep="")), nrow(coord_data_year))
  coord_data_year$end.date <- rep(as.Date(paste(year, "-06-30", sep="")), nrow(coord_data_year))
  coord_date_data <- rbind(coord_date_data, coord_data_year)
}

if (length(list.files("./data/modisdata/100kmsq/")) > 0){
  unaquired_coord_data <- UpdateSubsets(LoadDat = coord_date_data, StartDate = TRUE,
                                        Dir = "./data/modisdata/100kmsq/")
} else {
  unaquired_coord_data <- coord_date_data
}

MODISSubsets(LoadDat = unaquired_coord_data, Products = "MOD13Q1",
             Bands = c("250m_16_days_NDVI", "250m_16_days_pixel_reliability"),
             Size = c(5,5), SaveDir = "./data/modisdata/100kmsq/", StartDate = TRUE)
