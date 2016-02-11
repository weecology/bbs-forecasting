# Download and install PRISM data into PostgreSQL and values for BBS routes

library(prism)
library(raster)
library(dplyr)
library(tidyr)
library(stringr)
options(prism.path = "./data/prismdata")

#Years of prism data to 
yearsToUse=1981:2014
#Load the DB. create = TRUE does not seem 

sqliteDBFile='./data/bbsforecasting.sqlite'
if(file.exists(sqliteDBFile)){
  database <- src_sqlite("./data/bbsforecasting.sqlite", create = FALSE)
} else {
  database <- src_sqlite("./data/bbsforecasting.sqlite", create = TRUE)
}


#######################################################
#Downloads the raw prism rasters into the folder specified above. 
#######################################################
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
############################################################
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
#########################################################
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
  #As long as all of present==1 in the toCheck DF, then all prism rasters are accounted for.
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
########################################################
get_prism_data=function(){
  #Query the prism data and check to see if it throws an error. 
  prism_bbs_data=try(collect(tbl(database, sql('SELECT * from prism_bbs_data'))))
  
  #If it works it should be of type list. If it doesn't then it should be of class 'try-error'
  if(typeof(prism_bbs_data)=='list'){
    #If no error, return the data. 
    return(as.data.frame(prism_bbs_data))
  } else if(class(prism_bbs_data)=='try-error') { 
    #If the query returned an error, assume it was because the table doesn't exist
    #and need's to be created/loaded from raw prism data. 
    
    #Load the bbs data locations and convert them to a spatial object.
    #Stop here if bbs data isn't available. Could also make this query the DB as well.
    bbs_data <- try(read.csv("data/bbs_djata.csv"))
    if(class(bbs_data)=='try-error'){stop("Can't load bbs_data.csv inside get_prism_data()")}
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
    
  } else {
    
    stop('Im stopping because the sqlite query for prism_bbs data returned neither a data frame nor an error in get_prism_data()')
  }
  
}


#################################################################
#From raw prism monthly values calculate all the bioclim variables.
#You should not call this directly to load bioclim vars. Instead call
#get_bioclim_data(), which will 1st try to load the data from the sdqlite
#db before processing it all from scratch. 
###################################################################
#First a quick helper function to calculate some of the bioclim variables,
#like "precip in coldest month"
maxMinCombo=function(vec1,vec2,max=TRUE){
  #Return the value in vec1 in the position where
  #vec2 is either highest or lowest. But 1st check for na 
  #values. 
  if(any(is.na(vec1)) | any(is.na(vec2))){
    return(NA)
  } else  if(max){
    return(vec1[which.max(vec2)])
  } else {
    return(vec1[which.min(vec2)])
  }
}

process_bioclim_data=function(){
  #Get the prism data. It's pulled from the sqlite DB or downloaded as needed.
  prism_bbs_data=get_prism_data()
  
  #Spread out the climate variables ppt, tmean, etc into columns
  prism_bbs_data = prism_bbs_data %>%
    spread(clim_var, value)
  
  #Process the quarter ones first.
  quarterInfo=data.frame(month=1:12, quarter=c(1,1,1,2,2,2,3,3,3,4,4,4))
  bioclimQuarterData= prism_bbs_data %>%
    left_join(quarterInfo, by='month') %>%
    group_by(site_id, year, quarter) %>%
    summarize(precip=sum(ppt), temp=mean(tmean)) %>%
    ungroup() %>%
    group_by(site_id,year) %>%
    summarize(bio8=maxMinCombo(temp, precip, max=TRUE),
              bio9=maxMinCombo(temp, precip, max=FALSE),
              bio10=max(temp),
              bio11=min(temp),
              bio16=max(precip),
              bio17=min(precip),
              bio18=maxMinCombo(precip, temp, max=TRUE),
              bio19=maxMinCombo(precip, temp, max=FALSE)) %>%
    ungroup()
  
  #Next the yearly ones, joining the quartely ones  back in at the end. 
  bioclimData=prism_bbs_data %>%
    group_by(site_id, year) %>%
    mutate(monthlyTempDiff=tmax-tmin) %>%
    summarize(bio1=mean(tmean),
              bio2=mean(monthlyTempDiff),
              bio4=sd(tmean)*100,
              bio5=maxMinCombo(tmax,tmean,max=TRUE),
              bio6=maxMinCombo(tmin,tmean,max=FALSE),
              bio12=sum(ppt),
              bio13=max(ppt),
              bio14=min(ppt),
              bio15=cv(ppt)) %>%
    ungroup() %>%
    mutate(bio7=bio5-bio6,
           bio3=(bio2/bio7)*100) %>%
    full_join(bioclimQuarterData, by=c('site_id','year'))
  
  return(bioclimData)
}


####################################################################
#Load the bioclim variables from sqlite. if they aren't available,
#load the prism data, process bioclim, load it in sqlite, then return 
#bioclim as requested.

#Possibly put a check here to make sure all years are accounted for?
####################################################################

get_bioclim_data=function(){
  #Query the bioclim data and check to see if it throws an error. 
  bioclim_bbs_data=try(collect(tbl(database, sql('SELECT * from bioclim_bbs_data'))))
  
  #If it works it should be of type list. If it doesn't then it should be of class 'try-error'
  if(typeof(bioclim_bbs_data)=='list'){
    #If no error, return the data. 
    return(as.data.frame(bioclim_bbs_data))
  } else if(class(bioclim_bbs_data)=='try-error') { 
    #If erros on access the database, assume it's because the data just isn't there.
    #So load it from scratch, store it for future access, and return the newly process bioclim data
    bioclim_bbs_data=process_bioclim_data()
    
    x=copy_to(database, bioclim_bbs_data, temporary = FALSE, 
              indexes = list(c('site_id','year')))
     
    return(bioclim_bbs_data)
    
  } else {
    
    stop('Im stopping because the sqlite query for bioclim_bbs data returned neither a data frame nor an error in get_bioclim_data()')
  }
}