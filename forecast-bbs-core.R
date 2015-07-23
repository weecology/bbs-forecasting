# library(forecast)
library(dplyr)
library(broom)
library(MODISTools)
library(DBI)

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
  
  con <- dbConnect(RPostgres::Postgres(), dbname = 'postgres')
  data_path <- paste('./data/', 'bbs', '_data.csv', sep="")
  if (file.exists(data_path)){
    return(read.csv(data_path))
  }
  else{
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
  contig_pos <- na.contiguous(match(full_ts,unique(df$year)))
  filter(df, year %in% full_ts[contig_pos])
}

get_filtered_ts <- function(df, min_ts_length){
  data_by_site <- group_by(df, site_id)
  contig_ts <- do(data_by_site, get_longest_contig_ts(.))
  contig_ts_by_site <- group_by(contig_ts, site_id)
  contig_ts_length <- summarize(contig_ts_by_site, n_years = n_distinct(year))
  long_ts <- filter(contig_ts_length, n_years >= min_ts_length)
  contig_ts_long <- semi_join(contig_ts, long_ts)  
}

get_ndvi_ts_data <- function(tsdata){
  modis_data_files <- list.files("./data/modisdata/", pattern = "MODIS_Data", full.names = TRUE)
  if (length(modis_data_files) == 0){
    tsdata_by_site_yr <- group_by(tsdata, site_id, lat, long, year)
    richness_ndvi <- summarize(tsdata_by_site_yr, richness = n_distinct(species_id))
    richness_ndvi$start.date <- paste(richness_ndvi$year, "-06-01", sep = "")
    richness_ndvi$end.date <- paste(richness_ndvi$year, "-06-30", sep = "")
    
    output <- capture.output(
                 MODISSummaries(LoadDat = richness_ndvi, Dir = "./data/modisdata/",
                                Product = "MOD13Q1", Bands = c("250m_16_days_NDVI"),
                                ValidRange = c(-2000,10000),
                                NoDataFill = -3000, ScaleFactor = 0.0001, StartDate = TRUE)
    )
  }
  modis_data_files <- list.files("./data/modisdata/", pattern = "MODIS_Data", full.names = TRUE)
  richness_ndvi <- read.csv(modis_data_files[1])
  colnames(richness_ndvi)[9] <- "ndvi"
  return(richness_ndvi)
}