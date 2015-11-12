library(forecast)
library(dplyr)
library(broom)
library(MODISTools)
library(DBI)
library(ecoretriever)

database_exists <- function(schema_name, con){
  # Check to see if a database exists
  res <- dbSendQuery(con, "SELECT schema_name FROM information_schema.schemata")
  schema <- dbFetch(res)
  dbClearResult((res))
  return(schema_name %in% schema$schema_name)
}

install_dataset <- function(dataset){
  # Install a dataset using the ecoretriever
  
  #ecoretriever currently not working in RStudio, issue filed
  ecoretriever::install(dataset, 'postgres', conn_file = 'postgres_conn_file.txt')
}

get_bbs_data <- function(){
  # Get the BBS data
  
  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read.csv(data_path))
  }
  else{
    con <- dbConnect(RPostgres::Postgres(), dbname = 'postgres')
    if (!database_exists('bbs', con)){
      install_dataset('bbs')
    }

    #FIXME: This query doesn't currently deal with poorly sampled species (e.g., nocturnal)
    bbs_query = "SELECT (counts.statenum * 1000) + counts.route AS site_id, routes.lati as lat,
                        routes.loni as lon, counts.year, counts.aou AS species_id, counts.speciestotal AS abundance
                          FROM bbs.counts JOIN bbs.weather
                           ON bbs.counts.statenum=bbs.weather.statenum
                           AND bbs.counts.route=bbs.weather.route
                           AND bbs.counts.rpid=bbs.weather.rpid
                           AND bbs.counts.year=bbs.weather.year
                         JOIN bbs.routes
                           ON bbs.counts.statenum=bbs.routes.statenum
                           AND bbs.counts.route=bbs.routes.route
                         WHERE bbs.weather.runtype=1 AND bbs.weather.rpid=101;"
    bbs_results <- dbSendQuery(con, bbs_query)
    bbs_data <- dbFetch(bbs_results)
    colnames(bbs_data)[3] <- "long"
    write.csv(bbs_data, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_data)
  }
}

get_longest_contig_ts <- function(df){
  full_ts <- seq(from=min(df$year), to=max(df$year), by=1)
  years <- unique(df$year)
  contig_pos <- na.contiguous(match(full_ts, years))
  filter(df, year %in% years[contig_pos])
}

get_filtered_ts <- function(df, min_ts_length){
  data_by_site <- group_by(df, site_id)
  contig_ts <- do(data_by_site, get_longest_contig_ts(.))
  contig_ts_by_site <- group_by(contig_ts, site_id)
  contig_ts_length <- summarize(contig_ts_by_site, n_years = n_distinct(year))
  long_ts <- filter(contig_ts_length, n_years >= min_ts_length)
  contig_ts_long <- semi_join(contig_ts, long_ts)  
}

#' Add zeros to population dynamics data
#'
#' For population dynamics data that lacks real zeros, add them
#'
#' @param df data.frame including site_id, year, species_id, and abundance columns
#'
#' @return Data frame with real zeros added
get_popdyn_data_w_zeros <- function(commdyn_data, start_year, stop_year) {
  sp_years <- expand.grid(site_id = unique(commdyn_data$site_id), year = start_year:stop_year, species_id = unique(commdyn_data$species_id))
  popdyn_data <- commdyn_data %>%
    do(left_join(sp_years, ., by=c("site_id", "year", "species_id")))
  popdyn_data[is.na(popdyn_data)] <- 0
  return(popdyn_data)
}
get_ndvi_ts_data <- function(tsdata, modis_data_location){
  modis_data_files <- list.files(modis_data_location, pattern = "MODIS_Data", full.names = TRUE)
  if (length(modis_data_files) == 0){
    tsdata_by_site_yr <- group_by(tsdata, site_id, lat, long, year)
    richness_ndvi <- summarize(tsdata_by_site_yr, richness = n_distinct(species_id))
    richness_ndvi$start.date <- paste(richness_ndvi$year, "-06-01", sep = "")
    richness_ndvi$end.date <- paste(richness_ndvi$year, "-06-30", sep = "")
    
    output <- capture.output(
                 MODISSummaries(LoadDat = richness_ndvi, Dir = modis_data_location,
                                Product = "MOD13Q1", Bands = c("250m_16_days_NDVI"),
                                ValidRange = c(-2000,10000),
                                NoDataFill = -3000, ScaleFactor = 0.0001, StartDate = TRUE)
    )
  }
  modis_data_files <- list.files(modis_data_location, pattern = "MODIS_Data", full.names = TRUE)
  richness_ndvi <- read.csv(modis_data_files[1])
  richness_ndvi <- transmute(richness_ndvi, site_id = site_id, lat = lat, long = long,
                             year = year, richness = richness, start.date = start.date,
                             end.date = end.date, SubsetID = SubsetID,
                             ndvi = rowMeans(richness_ndvi[9:ncol(richness_ndvi)], na.rm = TRUE))
  return(richness_ndvi)
}

#' Perform a suite of time-series only forecasts
#'
#' Performs time-series only forecasting using naive, average, and ARIMA models
#' for each group in a grouped data frame of time-series.
#'
#' @param grouped_tsdata data.frame is a grouped data frame where the group
#' indicates how to split the data frame for running time-series models on each
#' group
#' @param timecol str name of the column that includes the time variable
#' @param responsecol str name of the column that includes the response variable
#' @param lag numeric how many time-steps to forecast and also how many
#'   time-steps in the current time-series to ignore when fitting.
#'
#' @return data.frame
#'
#' TODO: separate time-steps to hold out from fitting from time-steps to
#'       forecast, since when not hindcasting these may be different
get_ts_forecasts <- function(grouped_tsdata, timecol, responsecol, lag = 1){
  do(grouped_tsdata,
     year = .[[timecol]][(length(.[[responsecol]]) - lag + 1):length(.[[responsecol]])],
     cast_naive = naive(.[[responsecol]][1:(length(.[[responsecol]]) - lag)], lag),
     cast_avg = meanf(.[[responsecol]][1:(length(.[[responsecol]]) - lag)], lag),
     cast_arima = forecast(auto.arima(.[[responsecol]][1:(length(.[[responsecol]]) - lag)], seasonal = FALSE), h = lag),
     test_set = (.[[responsecol]][(length(.[[responsecol]]) - lag + 1):length(.[[responsecol]])])
  )
}

get_error_measures <- function(obs, pred){
  error <- obs - pred
  percenterror <- error / obs * 100
  me <- mean(error, na.rm=TRUE)
  rmse <- sqrt(mean(error^2, na.rm=TRUE))
  mae <- mean(abs(error), na.rm=TRUE)
  mpe <-  mean(percenterror, na.rm=TRUE)
  mape <- mean(abs(percenterror), na.rm=TRUE)
  results <- data.frame(t(unlist(c(me, rmse, mae, mpe, mape))))
  colnames(results) <- c('ME', 'RMSE', 'MAE', 'MPE', 'MAPE')
  results
}