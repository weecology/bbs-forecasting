# Download elevation data and place in SQLite database

#' Get route locations
#'
#' @param projection string projection for route data
#'
#' @return a spatial data frame including site_id, long, and lat
#' @importFrom sp SpatialPointsDataFrame
#' @importFrom dplyr collect copy_to src_sqlite src_tbls tbl %>%
get_route_data <- function(projection){
  bbs_data <- get_bbs_data()
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
#' @importFrom raster crs getData
get_elev_data <- function(){
  if (db_engine(action='check', table_to_check = 'elev_data')){
    return(db_engine(action='read', sql_query='SELECT * from elev_data'))
  } else {
    elevation <- getData("alt", country="US")[[1]]
    elev_proj <- crs(elevation)
    route_locations <- get_route_data(elev_proj)
    route_locations$elevs <- raster::extract(elevation, route_locations,
                                     buffer=40000, fun=mean) #buffers are in meters
    route_locations <- as.data.frame(route_locations)
    elev_data <- dplyr::select(route_locations, site_id, elevs)
    db_engine(action='write', df=elev_data, new_table_name = 'elev_data')
    return(elev_data)
  }
}
