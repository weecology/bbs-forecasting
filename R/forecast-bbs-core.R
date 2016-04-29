#' @importFrom DBI dbClearResult dbFetch dbSendQuery dbConnect
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


#' Filter poorly sampled BBS species
#'
#' Removes waterbirds, shorebirds, owls, kingfishers, knightjars,
#' dippers. These species are poorly sampled due to their aquatic or
#' noctural nature. Also removes taxa that were either partially unidentified
#' (e.g. "sp.") or were considered hybrids (e.g. "A x B").
#'
#' @param df dataframe containing an species_id column
#'
#' @return dataframe, filtered version of initial dataframe
#' @importFrom dplyr "%>%" inner_join do rowwise select filter group_by ungroup full_join n_distinct semi_join left_join
filter_species <- function(df){
  species_table = get_species_data()

  is_unidentified = function(names) {
    grepl("/|unid\\.|sp\\.| or |hybrid| X | x ", names)
  }

  valid_taxa = species_table %>%
    filter(!is_unidentified(english_common_name)) %>%
    filter(!is_unidentified(spanish_common_name)) %>%
    filter(aou > 2880) %>%
    filter(aou < 3650 | aou > 3810) %>%
    filter(aou < 3900 | aou > 3910) %>%
    filter(aou < 4160 | aou > 4210) %>%
    filter(aou != 7010)

  filter(df, species_id %in% valid_taxa$aou)
  }

combine_subspecies = function(df){

  species_table = get_species_data()

  # Subspecies have two spaces separated by non-spaces
  subspecies_names = species_table %>%
    filter(species_table$aou %in% unique(df$species_id)) %>%
    magrittr::extract2("spanish_common_name") %>%
    grep(" [^ ]+ ", ., value = TRUE)

  subspecies_ids = species_table %>%
    filter(spanish_common_name %in% subspecies_names) %>%
    magrittr::extract2("aou")

  # Drop the third word of the subspecies name to get the species name,
  # then find the AOU code
  new_subspecies_ids = species_table %>%
    slice(match(gsub(" [^ ]+$", "", subspecies_names),
                species_table$spanish_common_name)) %>%
    magrittr::extract2("aou")

  # replace the full subspecies names with species-level names
  for (i in seq_along(subspecies_ids)) {
    df$species_id[df$species_id == subspecies_ids[i]] = new_subspecies_ids[i]
  }

  df %>%
    group_by(site_id, year, species_id, lat, long) %>%
    summarize(abundance = sum(abundance)) %>%
    ungroup()
}

get_species_data = function() {
  data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
  if (file.exists(data_path)) {
    return(read.csv(data_path))
  }else{
    con <- dbConnect(RPostgres::Postgres(), dbname = 'postgres')
    species_table = dbFetch(dbSendQuery(con, "SELECT * FROM bbs.species;"))
    write.csv(species_table, file = data_path, row.names = FALSE, quote = FALSE)
    return(species_table)
  }
}

#' @export
get_bbs_data <- function(start_yr, end_yr, min_num_yrs){
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

    bbs_query = "SELECT (counts.statenum * 1000) + counts.route AS site_id, routes.latitude as lat,
                        routes.longitude as lon, counts.year, counts.aou AS species_id, counts.speciestotal AS abundance
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
    bbs_data <- dbFetch(bbs_results) %>%
      filter_species() %>%
      filter(year >= start_yr, year <= end_yr) %>%
      group_by(site_id) %>%
      filter(min(year) == start_yr, max(year) == end_yr, length(unique(year)) >= min_num_yrs) %>%
      combine_subspecies()
    colnames(bbs_data)[3] <- "long"
    write.csv(bbs_data, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_data)
  }
}

#' Get combined environmental data
#'
#' Master function for acquiring all environmental in a single table
#' @export
get_env_data <- function(){
  bioclim_data <- get_bioclim_data()
  elev_data <- get_elev_data()
  ndvi_data_raw <- get_bbs_gimms_ndvi()

  ndvi_data_summer <- ndvi_data_raw %>%
    filter(!is.na(ndvi), month %in% c('may', 'jun', 'jul'), year > 1981) %>%
    group_by(site_id, year) %>%
    dplyr::summarize(ndvi_sum = mean(ndvi))
  ndvi_data_winter <- ndvi_data_raw %>%
    filter(!is.na(ndvi), month %in% c('dec', 'jan', 'feb'), year > 1981) %>%
      group_by(site_id, year) %>%
        dplyr::summarize(ndvi_win = mean(ndvi))
  ndvi_data_ann <- ndvi_data_raw %>%
    filter(!is.na(ndvi), year > 1981) %>%
    group_by(site_id, year) %>%
    dplyr::summarize(ndvi_ann = mean(ndvi))
  ndvi_data <- inner_join(ndvi_data_summer, ndvi_data_winter, by = c('site_id', 'year'))
  ndvi_data <- inner_join(ndvi_data, ndvi_data_ann, by = c('site_id', 'year'))

  env_data <- full_join(bioclim_data, ndvi_data, by = c('site_id', 'year'))
  env_data <- full_join(env_data, elev_data, by = c('site_id'))
}

#' Get BBS population time-series data with environmental variables
#'
#' Selects sites with data spanning 1982 through 2013 containing at least 25
#' samples during that period.
#'
#' Attaches associated environmental data
#'
#' @param start_yr num first year of time-series
#' @param end_yr num last year of time-series
#' @param min_num_yrs num minimum number of years of data between start_yr & end_yr
#'
#' @return dataframe with site_id, lat, long, year, species_id, and abundance
#' @export
get_pop_ts_env_data <- function(start_yr, end_yr, min_num_yrs){
  bbs_data <- get_bbs_data(start_yr, end_yr, min_num_yrs)
  pop_ts_env_data <- bbs_data %>%
    add_env_data() %>%
    filter_ts(start_yr, end_yr, min_num_yrs) %>%
    dplyr::select(-lat, -long)
}

#' Get BBS richness time-series data with environmental variables
#'
#' @param start_yr num first year of time-series
#' @param end_yr num last year of time-series
#' @param min_num_yrs num minimum number of years of data between start_yr & end_yr
#'
#' @return dataframe with site_id, year, and richness
#' @export
get_richness_ts_env_data <- function(start_yr, end_yr, min_num_yrs){
  bbs_data <- get_bbs_data(start_yr, end_yr, min_num_yrs)
  richness_data <- bbs_data %>%
    group_by(site_id, year) %>%
    dplyr::summarise(richness = n_distinct(species_id)) %>%
    ungroup() %>%
    complete(site_id, year)
  richness_ts_env_data <- richness_data %>%
    add_env_data() %>%
    filter_ts(start_yr, end_yr, min_num_yrs)
}

#' Add environmental data to BBS data frame
#'
#' @param bbs_data dataframe that contains BBS site_id and year columns
#'
#' @return dataframe with original data and associated environmental data
add_env_data <- function(bbs_data){
  env_data <- get_env_data()
  bbs_data_w_env <- bbs_data %>%
    inner_join(env_data, by = c("site_id", "year"), copy = TRUE)
}

#' Filter BBS to specified time series period and number of samples
#'
#' @param bbs_data dataframe that contains BBS site_id and year columns
#' @param start_yr num first year of time-series
#' @param end_yr num last year of time-series
#' @param min_num_yrs num minimum number of years of data between start_yr & end_yr
#'
#' @return dataframe with original data and associated environmental data
filter_ts <- function(bbs_data, start_yr, end_yr, min_num_yrs){
  filterd_data <- bbs_data %>%
    filter(!is.na(bio1), !is.na(ndvi_sum), !is.na(ndvi_win), !is.na(elevs)) %>%
    group_by(site_id) %>%
    filter(min(year) == start_yr, max(year) == end_yr, length(unique(year)) >= min_num_yrs) %>%
    ungroup()
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
  contig_ts_length <- dplyr::summarize(contig_ts_by_site, n_years = n_distinct(year))
  long_ts <- filter(contig_ts_length, n_years >= min_ts_length)
  contig_ts_long <- semi_join(contig_ts, long_ts)
}

#' Get data for all contiguous samples ending in a given year
#'
#' Get data for all contiguous samples ending in a given year.
#' If only one year of contiguous data are available keep that year.
#'
#' @param df data.frame including a year column
#' @param stopyear year of the most recent sample
#'
#' @return Contiguous time-series data frame (data.frame).
#' A version of the original data frame only including years that
#' are in a contiguous series with stopyear.
get_modern_contig_ts <- function(df, stopyear) {

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

#' Cleanup output from multi-time-series forecasting into a simple data frame
#'
#' Takes the output of a dplyr based run of forecasting on multiple groups and
#' makes it into a simple data frame for further analysis. Currently this is
#' specific to the form of data produce insdie of get_ts_forecasts.
#'
#' @param ts_forecast_df data.frame containing year, site, cast_naive, cast_avg,
#'   and cast_arima columns
#'
#' @return data.frame
#' @export
cleanup_multi_ts_forecasts <- function(ts_forecast_df, groups){
  cleaned <-
    ts_forecast_df %>%
    rowwise() %>%
    do(
      {timeperiod <- as.data.frame(.$timeperiod)
      testset <- as.data.frame(.$test_set)
      colnames(timeperiod) <- c('timeperiod')
      colnames(testset) <- c('obs')
      naive <- as.data.frame(.$cast_naive)
      colnames(naive) <- c('pt_fcast', 'lo80', 'hi80', 'lo95', 'hi95')
      naive$model <- 'naive'
      naive <- cbind(timeperiod, naive, testset)
      avg <- as.data.frame(.$cast_avg)
      colnames(avg) <- c('pt_fcast', 'lo80', 'hi80', 'lo95', 'hi95')
      avg$model <- 'avg'
      avg <-cbind(timeperiod, avg, testset)
      arima <- as.data.frame(.$cast_arima)
      colnames(arima) <- c('pt_fcast', 'lo80', 'hi80', 'lo95', 'hi95')
      arima$model <- 'arima'
      arima <- cbind(timeperiod, arima, testset)
      df <- rbind(naive, avg, arima)
      for (group in groups) {df[[group]] <- .[[group]]}
      df %>% select(model, timeperiod, obs, everything())
      }
    )
  return(cleaned)
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
#' @importFrom forecast auto.arima forecast meanf naive
#' @export
get_ts_forecasts <- function(grouped_tsdata, timecol, responsecol, exogcol, lag = 1){
  get_train_data <- function(df, colname) df[[colname]][1:(nrow(df) - lag)]
  get_test_data <- function(df, colname) df[[colname]][(nrow(df) - lag + 1):nrow(df)]
  do(grouped_tsdata,
     timeperiod = get_test_data(., timecol),
     cast_naive = naive(get_train_data(., responsecol), lag),
     cast_avg = meanf(get_train_data(., responsecol), lag),
     cast_arima = forecast(auto.arima(get_train_data(., responsecol), seasonal = FALSE), h = lag),
     cast_exog_arima = forecast(auto.arima(get_train_data(., responsecol), xreg = get_train_data(., exogcol), seasonal = FALSE), xreg = get_test_data(., exogcol)),
     test_set = get_test_data(., responsecol)
  )
}

#' @export
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
