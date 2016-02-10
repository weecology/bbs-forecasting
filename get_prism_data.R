# Download and install PRISM data into PostgreSQL and values for BBS routes

library(prism)
library(raster)
library(dplyr)
library(tidyr)
library(stringr)
options(prism.path = "./data/prismdata")

yearsToUse=1981:2014

#######################################################
#Downloads the raw prism rasters into the folder specified above. 
downloadPrism=function(){
  months <- c(1:12)
  clim_vars <- c("ppt", "tmin", "tmean", "tmax")
  for (month in months){
    for (clim_var in clim_vars){
      get_prism_monthlys(type=clim_var, year = yearsToUse, month = month, keepZip=F)
    }
  }
}


#########################################################
#Pretty postgres isn't used anywhere after adding the raw prism to it,
#so i'll just comment this out for now. 
#datadirs = dir(datapath)
#for (datadir in datadirs) {
#  bil_file = paste(datadir, '.bil', sep = "")
#  bil_file_path = file.path(datapath, datadir, bil_file)
#  sql_file = paste(datadir, ".sql", sep = "")
#  sql_file_path = file.path(datapath, datadir, sql_file)
#  system(paste("raster2pgsql -s 4326", bil_file_path, datadir, ">", sql_file_path))
#  system(paste("psql -d bbsforecasting -f", sql_file_path))
#}

########################################################
#Takes output of ls_prism_data() and makes sure all files
#within a year range are there. 
check_If_Prism_Files_Present=function(prismLS, years){
  #Extract out all the variable names and dates
  prismLS = prismLS %>%
    rowwise() %>%
    mutate(var=strsplit(files, '_')[[1]][2], yearMonth=strsplit(files, '_')[[1]][5]) %>%
    mutate(year=as.integer(substr(yearMonth,1,4)), month=as.integer(substr(yearMonth, 5,6)), present=1) %>%
    select(-files, -yearMonth)
  
  #Setup a list of what should be there
  toCheck=expand.grid(year=years, month=1:12, var=c('ppt','tmax','tmean','tmin'))
  
  #Left join on the toCheck DF puts present=1 wherever that var/year/month combo was in the prismLS DF.
  #As long as all of present==1 in toCheck, then all prism rasters are accounted for.
  toCheck=toCheck %>%
    left_join(prismLS, by=c('year','month','var'))
  
  if(sum(toCheck$present, na.rm=TRUE)==nrow(toCheck)){
    return(TRUE)
  } else {
    
    return(FALSE)
  }
}

########################################################
#Return prism data for each site from the sqlite DB. If it's not loaded in
#the DB, then extract values from the raw prism rasters, store in DB, and return them.

get_prism_data=function(){
  database <- src_sqlite("./data/bbsforecasting.sqlite", create = TRUE)
  #Query the prism data and check to see if it throws an error. 
  prism_bbs_data=try(tbl(database, sql('SELECT * from prism_bbs_data')))
  if(!class(prism_bbs_data)=='try-error'){
    #If no error, return the data. 
    return(prism_bbs_data)
  } else { 
    #If the query returned an error, assume it was because the table doesn't exist
    #and need's to be created/loaded from raw prism data. 
    
    #Load the bbs data locations and convert them to a spatial object.
    bbs_data <- read.csv("data/bbs_data.csv")
    locations <- unique(select(bbs_data, site_id, long, lat))
    coordinates(locations) <- c("long", "lat")
    
    #Check to see if all the raw data in the years specified are downloaded,
    #download everything again if not. (Getting only what is needed, say if a 
    #previous download failed, might be overly complicated)
    if(!check_If_Prism_Files_Present(ls_prism_data(), yearsToUse)){
      downloadPrism()
    }
    
    #Load the prism data and extract using the bbs locations.
    prism_stacked <- prism_stack(ls_prism_data())
    extracted <- raster::extract(prism_stacked, locations)
    prism_bbs_data <- data.frame(site_id = locations$site_id, coordinates(locations), extracted)
    prism_bbs_data <- prism_bbs_data %>%
      gather(date, value, 4:ncol(prism_bbs_data)) %>%
      tidyr::extract(date, c("clim_var", "year", "month"),
                     "PRISM_([:alpha:]*)_stable_[:alnum:]*_([:digit:]{4})([:digit:]{2})_")
    
    #Format the data a little and load into the sqlite database.
    prism_bbs_data$year <- as.numeric(prism_bbs_data$year)
    prism_bbs_data$month <- as.numeric(prism_bbs_data$month)
    mydata <- copy_to(database, prism_bbs_data, temporary = FALSE,
                      indexes = list(c("site_id", "year", "month")))
    
    #Now return the data as asked for
    return(prism_bbs_data)
    
  }
  
}

