# library(forecast)
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
  
  con <- dbConnect(RPostgres::Postgres(), dbname = 'postgres')
  data_path <- paste('./data/', 'bbs', '_data.csv')
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
    write.csv(file = data_path, x = bbs_data)
    return(bbs_data)
  }
}
