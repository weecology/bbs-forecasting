#Functions for retrieving and processing NDVI data for BBS route locations. Currently only GIMMS is supported, with years 1981-2013 available.
#
#GIMMS: Call get_bbs_gimms_ndvi() to retrieve a dataframe of (site_id, year, month, ndvi). 
#-It will take a while to download the 14gb of data and extract values the 1st time around. Results are stored in an sqlite db for quick access.
#-To redo it just delete the database file.  
#-na values are due to missing periods from filtering out unwanted quality flags. 

library(rgdal)
library(gimms)
library(dplyr)
library(tidyr)
library(sp)
library(raster)
#library(doParallel)

sqlite_db_file='./data/bbsforecasting.sqlite'
database <- src_sqlite(sqlite_db_file, create = TRUE)
gimms_folder='./data/gimms_ndvi/'
dir.create(gimms_folder, showWarnings = FALSE, recursive = TRUE)
years_to_use=1981:2014

#################################################
#This assumes we want all gimms files that are available. It queries for files
#that are available for download and compares against files already in the gimms
#data directory. 
#Returns a list of files to download, which may be length 0. 
##################################################
get_gimms_download_list=function(){
  available_files_download_path=gimms::updateInventory()
  available_files_name=basename(available_files_download_path)
  
  files_present=list.files(gimms_folder)
  #hdr files are created from some of the gimms processing that we don't want to
  #use here.
  files_present=files_present[!grepl('hdr', files_present)]
  
  to_download=available_files_download_path[! available_files_name %in% files_present]
  
  return(to_download)
}

################################################
#Extract values from a single gimms file given a set of coordinates
################################################
extract_gimms_data=function(gimms_file_path, route_locations){
  #Have to load and extract twice. Once for the actual NDVI, once for the quality flag.
  gimmsRaster=rasterizeGimms(gimms_file_path, flag=FALSE)
  ndvi=raster::extract(gimmsRaster, route_locations)
  gimmsRaster=rasterizeGimms(gimms_file_path, flag=TRUE)
  flag=raster::extract(gimmsRaster, route_locations)
  
  year=as.numeric(substr(basename(gimms_file_path), 4,5))
  month=substr(basename(gimms_file_path), 6,8)
  day=substr(basename(gimms_file_path), 11,11)
  
  #Convert the a b to the 1st and 15th
  day=ifelse(day=='a',1,15)
  
  #Convert 2 digit year to 4 digit year
  year=ifelse(year>50, year+1900, year+2000)
  
  return(data.frame(year=year, month=month, day=day, ndvi=ndvi, flag=flag, site_id=route_locations@data$site_id, stringsAsFactors = FALSE))
}

################################################
#Extract the NDVI time series for all bbs routes
#from all years of gimms data
################################################
process_gimms_ndvi_bbs=function(){

  bbs_data <- try(read.csv("data/bbs_data.csv"))
  if(class(bbs_data)=='try-error'){stop("Can't load bbs_data.csv inside process_gimms_ndvi_bbs()")}
  route_locations <- unique(dplyr::select(bbs_data, site_id, long, lat))
  coordinates(route_locations) <- c("long", "lat")
  
  gimms_files=list.files(gimms_folder, full.names = TRUE)
  #hdr files are created from some of the gimms processing that we don't want to
  #use here.
  gimms_files=gimms_files[!grepl('hdr', gimms_files)]
  
  gimms_ndvi_bbs=data.frame()
  for(file_path in gimms_files){
    gimms_ndvi_bbs=extract_gimms_data(file_path, route_locations) %>%
      bind_rows(gimms_ndvi_bbs)
  }
  
  return(gimms_ndvi_bbs)
}

#################################################
#Filter the quality flags that we don't want.
#Average the bimonthly values into 1 value per month
#and add in -99 for na
#From the GIMMS readme:
#FLAG = 7 (missing data)
#FLAG = 6 (NDVI retrieved from average seasonal profile, possibly snow)
#FLAG = 5 (NDVI retrieved from average seasonal profile)
#FLAG = 4 (NDVI retrieved from spline interpolation, possibly snow)
#FLAG = 3 (NDVI retrieved from spline interpolation)
#FLAG = 2 (Good value)
#FLAG = 1 (Good value)
#################################################
filter_gimms_data=function(df){
  df=df %>%
    filter(flag<=3) %>%
    group_by(site_id, year, month) %>%
    summarize(ndvi=mean(ndvi)) %>%
    ungroup() %>%
    right_join( expand(df, site_id, month, year)) 
  return(df)
}


#################################################
#Get the GIMMS AVHRR ndvi bi-monthly time series for every bbs site.
#Pulling from the sqlite DB or extracting it from raw gimms data if needed. 
#################################################
get_bbs_gimms_ndvi = function(){
  if('gimms_ndvi_bbs_data' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from gimms_ndvi_bbs_data')), n = Inf))
  } else {
    print('Gimms NDVI bbs data not found, processing from scratch')
    
    files_to_download=get_gimms_download_list()
    if(length(files_to_download)>0){
      print('Downloading GIMMS data')
      downloadGimms(x=files_to_download, dsn=gimms_folder)
    }
    
    gimms_ndvi_bbs_data=process_gimms_ndvi_bbs()
    
    gimms_ndvi_bbs_data=filter_gimms_data(gimms_ndvi_bbs_data)
    
    copy_to(database, gimms_ndvi_bbs_data, temporary = FALSE, 
            indexes = list(c('site_id','year','month')))
    
    return(gimms_ndvi_bbs_data)
  }
}
