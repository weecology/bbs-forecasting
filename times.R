library(rgdal)
library(sp)
library(dplyr)
source("forecast-bbs-core.R")

locations = get_bbs_data() %>%
  dplyr::select(site_id, lat, long, year) %>%
  distinct()

#' @importFrom rgdal readOGR
#' @importFrom sp over SpatialPoints
get_time_zones = function() {
  if(!file.exists("data/tz_us.zip")) {
    download.file("http://efele.net/maps/tz/us/tz_us.zip", "data/tz_us.zip")
    unzip("data/tz_us.zip", exdir = "data/tz_us")
  }
  readOGR("data/tz_us/us/" ,"tz_us", verbose = FALSE)
}

add_time_zones = function(locations){
  boundaries = get_tz()
  time_zone = over(
    SpatialPoints(locations[, c("long", "lat")],
                  proj4string = CRS(proj4string(boundaries))),
    boundaries
  )
  cbind(locations, time_zone = time_zone[[1]])
}
