#library(rgdal)
library(gimms)
library(dplyr)
library(sp)
library(raster)
#library(doParallel)

sqlite_db_file='./data/bbsforecasting.sqlite'
database <- src_sqlite(sqlite_db_file, create = TRUE)
gimms_folder='./data/gimms_ndvi/'
years_to_use=1981:2014

#################################################
#This assumes we want all gimms files that are available. It queries for files
#that are available for download and checks to see which, if any, are missing.
#Returns TRUE if all files that are available are already present, and list of 
#files to download if any are missing.
##################################################
check_if_gimms_files_present=function(){
  available_files=sort(gimms::updateInventory()[1:3])
  available_files=basename(available_files)
  
  files_present=sort(list.files(gimms_folder))
  #hdr files are created from some of the gimms processing that we don't want to
  #use here.
  files_present=files_present[!grepl('hdr', files_present)]
  
  to_download=available_files[! available_files == files_present]
  
  if(length(to_download)==0){
    return(TRUE)
  } else {
    return(to_download)
  }
}


################################################
#Extract values from a single gimms file given a set of coordinates
################################################
extract_gimms_data=function(gimms_file_path, route_locations){
  #Have to load and extract twice. Once for the actual NDVI, once for the quality flag.
  gimmsRaster=rasterizeGimms(gimms_file_path, flag=FALSE)
  ndvi=extract(gimmsRaster, route_locations)
  gimmsRaster=rasterizeGimms(gimms_file_path, flag=TRUE)
  flag=extract(gimmsRaster, route_locations)
  
  year=as.numeric(substr(basename(gimms_file_path), 4,5))
  month=substr(basename(gimms_file_path), 6,8)
  day=substr(basename(gimms_file_path), 11,11)
  
  #Convert the a b to the 1st and 15th
  day=ifelse(day=='a',1,15)
  
  #Convert 2 digit year to 4 digit year
  year=ifelse(year>50, year+1900, year+2000)
  
  return(data.frame(year=year, month=month, day=day, ndvi=ndvi, flag=flag, site_id=route_locations@data$site_id))
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
#Get the GIMMS AVHRR ndvi bi-monthly time series for every bbs site.
#Pulling from the sqlite DB or extracting it from raw gimms data
#################################################
get_bbs_gimms_ndvi = function(){
  if('gimms_ndvi_bbs_data' %in% src_tbls(database)){
    return(collect(tbl(database, sql('SELECT * from gimms_ndvi_bbs_data'))))
  } else {
    print('Gimms NDVI bbs data not found, processing from scratch')
    
    file_status=check_if_gimms_files_present(gimms_folder)
    if(!isTRUE(file_status)){
      print('Downloading GIMMS data')
      downloadGimms(x=file_status, dsn=gimms_folder)
    }
    
    gimms_ndvi_bbs_data=process_gimms_ndvi_bbs()
    
    copy_to(database, gimms_ndvi_bbs_data, temporary = FALSE, 
            indexes = list(c('site_id','year','month','day')))
    
    return(bioclim_bbs_data)
    
  }
}










# Get NDVI data from MODIS

library(MODISTools)
library(dplyr)

route_data <- read.csv('data/BBS_routes.csv')
route_data$site_id <-  route_data$statenum * 1000 + route_data$route

contig_modern_routes <- read.csv('data/contig_modern_routes.csv')

colnames(contig_modern_routes) <- c('site_id')
route_data <- semi_join(route_data, contig_modern_routes)


coord_data <- route_data[c('lati', 'loni')]
coord_data <- unique(coord_data) # the Routes table has duplicate coordinate values for some reason
colnames(coord_data) <- c('lat', 'long')
coord_date_data <- coord_data
coord_date_data$start.date <- rep(as.Date("2000-06-01"), nrow(coord_data))
coord_date_data$end.date <- rep(as.Date("2000-06-30"), nrow(coord_data))

years = 2001:2015
for (year in years){
  coord_data_year <- coord_data
  coord_data_year$start.date <- rep(as.Date(paste(year, "-06-01", sep="")), nrow(coord_data_year))
  coord_data_year$end.date <- rep(as.Date(paste(year, "-06-30", sep="")), nrow(coord_data_year))
  coord_date_data <- rbind(coord_date_data, coord_data_year)
}

if (length(list.files("./data/modisdata/100kmsq/")) > 0){
  unaquired_coord_data <- UpdateSubsets(LoadDat = coord_date_data, StartDate = TRUE,
                                        Dir = "./data/modisdata/100kmsq/")
} else {
  unaquired_coord_data <- coord_date_data
}

MODISSubsets(LoadDat = unaquired_coord_data, Products = "MOD13Q1",
             Bands = c("250m_16_days_NDVI", "250m_16_days_pixel_reliability"),
             Size = c(5,5), SaveDir = "./data/modisdata/100kmsq/", StartDate = TRUE)
