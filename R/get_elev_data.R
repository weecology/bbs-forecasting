# Download elevation data and place in SQLite database

#' Get Elevation data
#'
#' If elevation data already exists, do nothing, otherwise use the raster
#' package to download the average elevation in a 40km radius of each BBS route
#' and push to the SQLite database
#'
#' @return a data frame including site_id and elevs columns
#' @importFrom raster getData
get_elev_data <- function(){
  if (db_engine(action='check', table_to_check = 'elev_data')){
    return(db_engine(action='read', sql_query='SELECT * from elev_data'))
  } else {
    dir.create("data/alt", showWarnings = FALSE)
    elevation <- getData("alt", country="US", path="data/alt")[[1]]
    route_locations <- get_route_data()
    route_locations$elevs <- raster::extract(elevation, route_locations,
                                     buffer=40000, fun=mean) #buffers are in meters
    route_locations <- as.data.frame(route_locations)
    elev_data <- dplyr::select(route_locations, site_id, elevs)
    db_engine(action='write', df=elev_data, new_table_name = 'elev_data')
    return(elev_data)
  }
}
