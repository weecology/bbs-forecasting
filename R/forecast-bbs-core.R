#' Install a particular dataset via rdataretreiver
#'
#' @param dataset name
install_dataset <- function(dataset){
  # Install a dataset using the rdataretriever

  rdataretriever::install(dataset, 'sqlite', db_file='./data/bbsforecasting.sqlite')
}

#' Single wrapper for all database actions
#'
#' We require only a few simple sql methods. They are 1. Writing an entire dataframe
#' directly to a database as it's own table, 2. Reading the same tables as dataframes,
#' possibly with some modification using SQL statements, 3. Checking to see if a
#' table exists. If a particular table does it exists it's assumed it has all data
#' required.
#'
#' read returns a dataframe
#' write returns nothing
#' check returns boolean
#'
#' @param action Action to perform in db call. Either read, write, or check
#' @param db name of database. A file if using sqlite
#' @param sql_query SQL statement if action is read
#' @param df Dataframe of data if action is write. Will copy the dataframe verbatim to it's own table with name new_table_name
#' @param new_table_name Table name for new data being written
#' @param table_to_check Table name to check if it exists for when action is check
#' @importFrom dplyr copy_to src_sqlite src_tbls collect tbl

db_engine=function(action, db='./data/bbsforecasting.sqlite', sql_query=NULL,
                   df=NULL, new_table_name=NULL, table_to_check=NULL){

  if(!dir.exists("data")){dir.create("data")}

  con <- src_sqlite(db, create=TRUE)

  if(action=='read'){
    to_return=collect(tbl(con, sql(sql_query)), n=Inf)

  } else if(action=='write') {
    copy_to(con, df, name=new_table_name, temporary = FALSE)
    to_return=NA

  } else if(action=='check') {
    #Only works with sqlite for now.
    to_return=tolower(table_to_check) %in% tolower(src_tbls(con))

  } else {
    stop(paste0('DB action: ',action,' not found'))
  }

  #Close the connection before returning results.
  rm(con)
  return(to_return)
}

#' Filter poorly sampled BBS species
#'
#' Removes waterbirds, shorebirds, owls, kingfishers, knightjars,
#' dippers. These species are poorly sampled due to their aquatic or
#' noctural nature. Also removes taxa that were either partially unidentified
#' (e.g. "sp.") or were considered hybrids (e.g. "A x B") or were listed as more
#' than one species (e.g. "A / B")
#'
#' @param df dataframe containing an species_id column
#'
#' @return dataframe, filtered version of initial dataframe
#' @importFrom dplyr "%>%" inner_join do rowwise select filter group_by ungroup full_join n_distinct semi_join left_join
filter_species <- function(df){
  species_table = get_species_data()

  is_unidentified = function(names) {
    #Befor filter account for this one hybrid of 2 subspecies so it's kept
    names[names=='auratus auratus x auratus cafer']='auratus auratus'
    grepl('sp\\.| x |\\/', names)
  }

  valid_taxa = species_table %>%
    filter(!is_unidentified(species)) %>%
    filter(AOU > 2880) %>%
    filter(AOU < 3650 | AOU > 3810) %>%
    filter(AOU < 3900 | AOU > 3910) %>%
    filter(AOU < 4160 | AOU > 4210) %>%
    filter(AOU != 7010)

  filter(df, species_id %in% valid_taxa$AOU)
}

#' Combine subspecies into their common species
#'
#' @importFrom dplyr "%>%" filter slice group_by summarise ungroup
#' @importFrom magrittr extract2
#' @importFrom stringr word
combine_subspecies = function(df){

  species_table = get_species_data()

  # Subspecies have two spaces separated by non-spaces
  subspecies_names = species_table %>%
    filter(species_table$AOU %in% unique(df$species_id)) %>%
    magrittr::extract2("spanish_common_name") %>%
    grep(" [^ ]+ ", ., value = TRUE)

  subspecies_ids = species_table %>%
    filter(spanish_common_name %in% subspecies_names) %>%
    extract2("AOU")

  # Drop all but the first two words to get the root species name,
  # then find the AOU code
  new_subspecies_ids = species_table %>%
    slice(match(word(subspecies_names, 1,2),
                species_table$spanish_common_name)) %>%
    extract2("AOU")

  # replace the full subspecies names with species-level names
  for (i in seq_along(subspecies_ids)) {
    df$species_id[df$species_id == subspecies_ids[i]] = new_subspecies_ids[i]
  }

  df %>%
    group_by(site_id, year, species_id, lat, long) %>%
    summarise(abundance = sum(abundance)) %>%
    ungroup()
}

get_species_data = function() {
  data_path <- paste('./data/', 'bbs', '_species.csv', sep = "")
  if (file.exists(data_path)) {
    return(read.csv(data_path))
  }else{
    species_table=db_engine(action = 'read', sql_query = 'SELECT * FROM breed_bird_survey_species')
    write.csv(species_table, file = data_path, row.names = FALSE, quote = FALSE)
    save_provenance(species_table)
    return(species_table)
  }
}

#' Get the primary bbs data file which compiles the counts, route info, and weather
#' data. Install it via rdataretriever if needed.
#'
#' @export
#' @importFrom dplyr "%>%" group_by
#' @importFrom readr read_csv
get_bbs_data <- function(){

  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read_csv(data_path))
  }
  else{

    if (!db_engine(action='check', table_to_check = 'breed_bird_survey_counts')){
      install_dataset('breed-bird-survey')
    }

    #Primary BBS dataframe
    bbs_query ="SELECT
                  (counts.statenum*1000) + counts.Route AS site_id,
                  Latitude AS lat,
                  Longitude AS long,
                  AOU AS species_id,
                  counts.Year AS year,
                  speciestotal AS abundance
                FROM
                  breed_bird_survey_counts AS counts
                  JOIN breed_bird_survey_weather
                    ON counts.statenum=breed_bird_survey_weather.statenum
                    AND counts.route=breed_bird_survey_weather.route
                    AND counts.rpid=breed_bird_survey_weather.rpid
                    AND counts.year=breed_bird_survey_weather.year
                  JOIN breed_bird_survey_routes
                    ON counts.statenum=breed_bird_survey_routes.statenum
                    AND counts.route=breed_bird_survey_routes.route
                WHERE breed_bird_survey_weather.runtype=1 AND breed_bird_survey_weather.rpid=101"

    bbs_data=db_engine(action='read', sql_query = bbs_query) %>%
      filter_species() %>%
      group_by(site_id) %>%
      combine_subspecies()
    save_provenance(bbs_data)
    write.csv(bbs_data, file = data_path, row.names = FALSE, quote = FALSE)
    return(bbs_data)
  }
}

#' Get route locations in a SpatialPointsDataFrame
#'
#' @param projection string projection for route data
#'
#' @return a spatial data frame including site_id, long, and lat
#' @importFrom sp SpatialPointsDataFrame CRS
#' @importFrom dplyr collect copy_to src_sqlite src_tbls tbl %>%
get_route_data <- function(){
  p=CRS('+proj=longlat +datum=WGS84')
  bbs_data <- get_bbs_data()
  route_locations <- unique(dplyr::select(bbs_data, site_id, long, lat))
  spatial_routes <- route_locations %>%
    dplyr::select(long, lat) %>%
    SpatialPointsDataFrame(data=route_locations, proj4string=p)
}

#' Get combined environmental data
#'
#' Master function for acquiring all environmental in a single table
#' @export
#' @importFrom dplyr "%>%" filter group_by summarise ungroup inner_join full_join

get_env_data <- function(){
  bioclim_data <- get_bioclim_data()
  elev_data <- get_elev_data()
  ndvi_data_raw <- get_bbs_gimms_ndvi()
  
  #Offset the NDVI year by 6 months so that the window for will be July 1 - June 30. 
  #See https://github.com/weecology/bbs-forecasting/issues/114
  ndvi_data_raw$year = with(ndvi_data_raw, ifelse(month %in% c('jul','aug','sep','oct','nov','dec'), year+1, year))
  
  ndvi_data_summer <- ndvi_data_raw %>%
    filter(!is.na(ndvi), month %in% c('may', 'jun', 'jul'), year > 1981) %>%
    group_by(site_id, year) %>%
    summarise(ndvi_sum = mean(ndvi)) %>%
    ungroup()
  ndvi_data_winter <- ndvi_data_raw %>%
    filter(!is.na(ndvi), month %in% c('dec', 'jan', 'feb'), year > 1981) %>%
    group_by(site_id, year) %>%
    summarise(ndvi_win = mean(ndvi)) %>%
    ungroup()
  ndvi_data_ann <- ndvi_data_raw %>%
    filter(!is.na(ndvi), year > 1981) %>%
    group_by(site_id, year) %>%
    summarise(ndvi_ann = mean(ndvi)) %>%
    ungroup()
  ndvi_data <- inner_join(ndvi_data_summer, ndvi_data_winter, by = c('site_id', 'year'))
  ndvi_data <- inner_join(ndvi_data, ndvi_data_ann, by = c('site_id', 'year'))

  env_data <- full_join(bioclim_data, ndvi_data, by = c('site_id', 'year'))
  env_data <- full_join(env_data, elev_data, by = c('site_id'))
  save_provenance(env_data)
  return(env_data)
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
#' @importFrom dplyr "%>%" filter select
get_pop_ts_env_data <- function(start_yr, end_yr, min_num_yrs){
  pop_ts_env_data = get_bbs_data() %>%
    filter_ts(start_yr, end_yr, min_num_yrs) %>%
    complete(site_id, year) %>% 
    add_env_data() %>%
    filter(!is.na(bio1), !is.na(ndvi_sum), !is.na(elevs)) %>%
    group_by(site_id) %>%
    #Filter min_num_years again after accounting for missing environmental data
    filter(length(unique(year)) >= min_num_yrs) %>%
    ungroup() %>%
    add_observers()

  save_provenance(pop_ts_env_data)
  return(pop_ts_env_data)
}

#' Get BBS richness time-series data with environmental variables. Will
#' fill in NA values for richness in missing years as long as a site
#' meets the year requirments and has environmental data.
#'
#' @param start_yr num first year of time-series
#' @param end_yr num last year of time-series
#' @param min_num_yrs num minimum number of years of data between start_yr & end_yr
#'
#' @return dataframe with site_id, year, and richness
#' @export
#' @importFrom dplyr "%>%" left_join select distinct group_by summarise ungroup filter
#' @importFrom tidyr complete
get_richness_ts_env_data <- function(start_yr, end_yr, min_num_yrs){
  richness_ts_env_data = get_pop_ts_env_data(start_yr, end_yr, min_num_yrs) %>%
    collapse_to_richness()
  save_provenance(richness_ts_env_data)
  return(richness_ts_env_data)
}

#' Replace species-level information with richness values, eliminating redundant rows
#' @param df a data frame, such as produced by get_bbs_data or get_pop_ts_env_data
#' @export
#' @importFrom dplyr n group_by mutate select ungroup distinct right_join
collapse_to_richness = function(df){

  # Code below assumes that NA abundance always means no observations for the
  # whole transect
  stopifnot(all(is.na(df$abundance) == is.na(df$species_id)))

  df %>%
    group_by(site_id, year) %>%
    mutate(richness = sum(abundance > 0)) %>%
    select(-species_id, -abundance) %>% 
    ungroup() %>%
    distinct() %>% 
    right_join(distinct(select(df, -species_id, -abundance)))
}


#' Add environmental data to BBS data frame
#'
#' @param bbs_data dataframe that contains BBS site_id and year columns
#'
#' @return dataframe with original data and associated environmental data
#' @importFrom dplyr "%>%" inner_join
add_env_data <- function(bbs_data){
  env_data <- get_env_data()
  bbs_data_w_env <- bbs_data %>%
    inner_join(env_data, by = c("site_id", "year"), copy = TRUE)
}


#' Add observer ID's to data
#' @param bbs_data dataframe that contains BBS site_id and year columns
#' @importFrom dplyr %>% left_join
add_observers = function(bbs_data) {
  data_path = "data/bbs_observers.csv"
  if (file.exists(data_path)) {
    observer_info = read_csv(data_path)
  } else {
    observer_query = '
    SELECT
                    breed_bird_survey_weather.obsn as observer_id,
                    year,
                    (statenum*1000)+route as site_id
                  FROM
                    breed_bird_survey_weather
                  WHERE
                    runtype=1 AND rpid=101'

    observer_info = db_engine(action = 'read', sql_query = observer_query)
    write.csv(observer_info, file = data_path, row.names = FALSE, quote = FALSE)
  }  
  bbs_data %>%
    left_join(observer_info, by = c('site_id','year'))
}



#' Filter BBS to specified time series period and number of samples
#'
#' @param bbs_data dataframe that contains BBS site_id and year columns
#' @param start_yr num first year of time-series
#' @param end_yr num last year of time-series
#' @param min_num_yrs num minimum number of years of data between start_yr & end_yr
#'
#' @return dataframe with original data and associated environmental data
#' @importFrom dplyr "%>%" filter group_by summarise ungroup
filter_ts <- function(bbs_data, start_yr, end_yr, min_num_yrs){
  sites_to_keep = bbs_data %>%
    filter(year >= start_yr, year <= end_yr) %>%
    group_by(site_id) %>%
    summarise(num_years=length(unique(year))) %>%
    ungroup() %>%
    filter(num_years >= min_num_yrs)

  filterd_data <- bbs_data %>%
    filter(year >= start_yr, year <= end_yr) %>%
    filter(site_id %in% sites_to_keep$site_id)
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
get_ts_forecasts <- function(grouped_tsdata, timecol, responsecol, exogcols,
                             lag = 1, pred_int_levels = c(80, 95)){
  get_train_data <- function(df, colnames) as.matrix(df[, colnames])[1:(nrow(df) - lag),]
  get_test_data <- function(df, colnames) as.matrix(df[, colnames])[(nrow(df) - lag + 1):nrow(df),]
  tsmodel_forecasts <- do(grouped_tsdata,
     timeperiod = get_test_data(., timecol),
     cast_naive = naive(get_train_data(., responsecol), lag, level = pred_int_levels),
     cast_avg = meanf(get_train_data(., responsecol), lag, level = pred_int_levels),
     cast_arima = forecast(auto.arima(get_train_data(., responsecol), seasonal = FALSE), h = lag, level = pred_int_levels),
     cast_exog_arima = forecast(auto.arima(get_train_data(., responsecol), xreg = get_train_data(., exogcols), seasonal = FALSE), xreg = get_test_data(., exogcols), level = pred_int_levels),
     cast_ets = forecast(ets(get_train_data(., responsecol)), lag, level = pred_int_levels),
     test_set = get_test_data(., responsecol)
  )
  cleaned_ts_model_forecasts <- cleanup_ts_forecasts(tsmodel_forecasts)
  ts_fcasts <- restruct_ts_forecasts(cleaned_ts_model_forecasts)
  save_provenance(ts_fcasts)
  return(ts_fcasts)
}

#' Cleanup time-series forecast output
#'
#' Transforms a data frame with one row per for each site and columns
#' with cells containing full forecast output a number of different forecasts,
#' into one where each row contains the point forecast and forecast intervals
#' for a single site, time period, model combination.
#'
#' @param ts_model_forecasts dataframe with site_id, timeperiod, cast_naive,
#'   cast_arima, cast_exog_arima, and test_set columns. Each cast_ column
#'   contains the full forecast object for each year of the forecast
#'
#' @return data.frame
cleanup_ts_forecasts <- function(tsmodel_forecasts){
  column_names <-
  tsmodel_forecasts %>%
    do(
       {timeperiod <- as.data.frame(.$timeperiod)
        colnames(timeperiod) <- c('timeperiod')
        naive <- as.data.frame(.$cast_naive)
        naive <- cbind(timeperiod, naive)
        naive$model <- 'naive'
        naive$site_id <- .$site
        naive$obs <- .$test_set
        avg <- as.data.frame(.$cast_avg)
        avg <-cbind(timeperiod, avg)
        avg$model <- 'avg'
        avg$site_id <- .$site
        avg$obs <- .$test_set
        arima <- as.data.frame(.$cast_arima)
        arima <- cbind(timeperiod, arima)
        arima$model <- 'arima'
        arima$site_id <- .$site
        arima$obs <- .$test_set
        exog_arima <- as.data.frame(.$cast_exog_arima)
        exog_arima <- cbind(timeperiod, exog_arima)
        exog_arima$model <- 'exog_arima'
        exog_arima$site_id <- .$site
        exog_arima$obs <- .$test_set
	ets <- as.data.frame(.$cast_ets)
        ets <- cbind(timeperiod, ets)
        ets$model <- 'ets'
        ets$site_id <- .$site
        ets$obs <- .$test_set
        df <- rbind(naive, avg, arima, exog_arima, ets)
        df %>% dplyr::select(site_id, model, timeperiod, obs, everything())
       }
      )
}

#' Split time-series forecasts into point forecast and prediction interval dfs
#'
#' @param cleaned_ts_model_forecasts dataframe with site_id, model, timeperiod,
#'   obs, pt_fcast, and Lo and Hi columns for each level of prediction interval
#'
#' @return list of data.frames. The first value is the pt_fcast part of the
#'   input data frame. The second value is a long version of the full data.frame
#'   with columns for the level and interval characterizing each prediction
#'   interval
restruct_ts_forecasts <- function(cleaned_ts_model_fcasts){
  cleaned_ts_model_fcasts <- dplyr::rename(cleaned_ts_model_fcasts, pt_fcast = `Point Forecast`)
  ts_pt_fcasts <- dplyr::select(cleaned_ts_model_fcasts, site_id:pt_fcast)
  ts_fcasts_pred_int <- cleaned_ts_model_fcasts %>%
    tidyr::gather(level, interval_edge, -site_id, -model, -timeperiod, -obs, -pt_fcast) %>%
    tidyr::separate(level, into = c("hilo", "levels")) %>%
    dplyr::group_by(site_id, model, timeperiod, obs, pt_fcast, levels) %>%
    dplyr::summarize(lo = min(interval_edge), hi = max(interval_edge))
  ts_fcasts_pred_int$levels <- as.numeric(ts_fcasts_pred_int$levels)
  ts_fcasts <- list(pt_est = ts_pt_fcasts, intervals = ts_fcasts_pred_int)
  save_provenance(ts_fcasts)
  return(ts_fcasts)
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
  save_provenance(results)
  return(results)
}

#' Visualize a forecast
#'
#' Visualizes a forecast for a single focal site from a group of forecasts
#'
viz_forecast <- function(data, forecasts, focal_site){
  forecasts_focal <- filter(forecasts, site_id == focal_site)
  train_set <- filter(data, site_id == focal_site, year < min(forecasts$timeperiod))
  test_set <- filter(data, site_id == focal_site, year >= min(forecasts$timeperiod))

  examp_fcasts <- ggplot(train_set, aes(year, richness)) +
    geom_point(size = 2) +
    geom_line(size = 0.75) +
    geom_point(data = test_set, size = 2) +
    geom_line(data = test_set, linetype = "dashed", size = 0.75) +
    geom_line(data = forecasts_focal,
              aes(x = timeperiod, y = pt_fcast, color = model),
              size = 1) +
    scale_color_brewer(palette = "Set3") +
    labs(x = "Year", y = "Diversity")

  examp_fcasts
}
