# Download elevation data and place in SQLite database

library(raster)
library(sp)
library(dplyr)

#' Get route locations
#'
#' @param projection string projection for route data
#'
#' @return a spatial data frame including site_id, long, and lat
get_route_data <- function(projection){
  bbs_data <- try(read.csv("data/bbs_data.csv"))
  if(class(bbs_data)=='try-error'){stop("Can't load bbs_data.csv")}
  route_locations <- unique(dplyr::select(bbs_data, site_id, long, lat))
  spatial_routes <- route_locations %>%
    dplyr::select(long, lat) %>%
    SpatialPointsDataFrame(data=route_locations, proj4string=projection)
}

#' Get Elevation data
#'
#' If elevation data already exists, do nothing, otherwise use the raster
#' package to download the average elevation in a 40km radius of each BBS route
#' and push to the SQLite database
#'
#' @return a data frame including site_id and elevs columns
get_elev_data <- function(){
  sqlite_db_file='./data/bbsforecasting.sqlite'
  database <- src_sqlite(sqlite_db_file, create=TRUE)
  if('elev_data' %in% src_tbls(database)){
    return(data.frame(collect(tbl(database, "elev_data"))))
  } else { 
    elevation <- getData("alt", country="US")[[1]]
    elev_proj <- crs(elevation)
    route_locations <- get_route_data(elev_proj)
    route_locations$elevs <- raster::extract(elevation, route_locations,
                                     buffer=40000, fun=mean) #buffers are in meters
    route_locations <- as.data.frame(route_locations)
    elev_data <- dplyr::select(route_locations, site_id, elevs)
    copy_to(database, elev_data, temporary=FALSE, indexes = list(c("site_id")))
    return(elev_data)
  }
}

get_elev_data()
