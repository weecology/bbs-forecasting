library(rgdal)
library(sp)
library(dplyr)
library(lubridate)
library(maptools)
devtools::load_all()



#' @importFrom rgdal readOGR
#' @importFrom sp over SpatialPoints
get_time_zones = function() {
  if(!file.exists("data/tz_us.zip")) {
    download.file("http://efele.net/maps/tz/us/tz_us.zip", "data/tz_us.zip")
    unzip("data/tz_us.zip", exdir = "data/tz_us")
  }
  readOGR("data/tz_us/us/" ,"tz_us", verbose = FALSE)
}

add_time_zones = function(locations){
  boundaries = get_time_zones()
  time_zone = over(
    SpatialPoints(locations[, c("long", "lat")],
                  proj4string = CRS(proj4string(boundaries))),
    boundaries
  )
  out = cbind(locations, time_zone = time_zone[[1]])
  out$time_zone = as.character(out$time_zone)

  # Fix two coastal routes that were misclassified as NA
  out$time_zone[out$site_id %in% c(47007, 88029)] = "America/New_York"
  out
}

add_times = function(df){
  # BBS reports times as integers in HMM format.  So "855" is "08:55 AM"
  # Here, I record the time in the original time zone, then convert to UTC.
  # So the outputted times should be late morning or early afternoon in UTC time.
  df %>%
    mutate(date_time = paste(year, month, day, start_time %/% 100, start_time %% 100, sep = "-")) %>%
    mutate(date_time = ymd_hm(date_time, tz = time_zone[1])) %>%
    mutate(date_time = with_tz(date_time, tzone = "UTC"))
}

events = get_bbs_data() %>%
  dplyr::select(site_id, lat, long, year, month, day, start_time) %>%
  distinct() %>%
  add_time_zones() %>%
  group_by(time_zone) %>%
  do(add_times(.)) %>%
  ungroup() %>%
  mutate(yday = yday(date_time)) %>%
  dplyr::select(-start_time, -time_zone, -month, -day)

# From the crepuscule help:
# > Input can consist of one location and at least one POSIXct times, or one
# > POSIXct time and at least one location. solarDep is recycled as needed.
# Does that mean one of them needs to be a scalar??
# solarDep==6 is for "civil dawn", sun 6 degrees below horizon
x = crepuscule(
  cbind(events$long, events$lat),
  events$date_time,
  solarDep = 6,
  direction = "dawn",
  POSIXct.out = TRUE
) %>%
  bind_cols(events) %>%
  mutate(min_post_dawn = as.double(date_time - time, units = "mins")) %>%
  select(-day_frac, -time, -date_time)

# * Site 35039 is always off by almost exactly an hour 22 times & is right by an intra-state time zone change.
#
# * Site 35031 is also off by almot exactly an hour one time, and it's in Indiana, which
#   has super-complicated time zones & DST (https://en.wikipedia.org/wiki/Time_in_Indiana).
