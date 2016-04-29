#' @importFrom rgdal readOGR
get_time_zones = function() {
  if(!file.exists("data/tz_us.zip")) {
    download.file("http://efele.net/maps/tz/us/tz_us.zip", "data/tz_us.zip")
    unzip("data/tz_us.zip", exdir = "data/tz_us")
  }
  readOGR("data/tz_us/us/" ,"tz_us", verbose = FALSE)
}

#' @importFrom sp proj4string SpatialPoints CRS over
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

#' @importFrom lubridate yday day month year ymd_hm
add_times = function(df){
  # BBS reports times as integers in HMM format.  So "855" is "08:55 AM"
  # Here, I record the time in the original time zone, then convert to UTC.
  # So the outputted times should be late morning or early afternoon in UTC time.
  df %>%
    mutate(date_time = paste(year, month, day, start_time %/% 100, start_time %% 100, sep = "-")) %>%
    mutate(date_time = ymd_hm(date_time, tz = time_zone[1])) %>%
    mutate(date_time = with_tz(date_time, tzone = "UTC"))
}

#' @importFrom dplyr bind_cols do distinct
#' @importFrom maptools crepuscule
get_temporal_data = function(start_yr, end_yr, min_num_yrs){
  events = get_bbs_data(start_yr = start_yr, end_yr = end_yr, min_num_yrs = min_num_yrs) %>%
    select(site_id, lat, long, year, month, day, start_time) %>%
    distinct() %>%
    add_time_zones() %>%
    group_by(time_zone) %>%
    do(add_times(.)) %>%
    ungroup() %>%
    mutate(yday = yday(date_time)) %>%
    dplyr::select(-start_time, -time_zone, -month, -day)

  # solarDep==6 is for "civil dawn", sun 6 degrees below horizon
  out = crepuscule(
    cbind(events$long, events$lat),
    events$date_time,
    solarDep = 6,
    direction = "dawn",
    POSIXct.out = TRUE
  ) %>%
    bind_cols(events) %>%
    mutate(min_post_dawn = as.double(date_time - time, units = "mins")) %>%
    select(-day_frac, -time, -date_time)

  # Site 35039 is off by almost exactly an hour 22 times & is right by an
  # intra-state time zone change.
  is_hour_early = out$site_id == 35039 & out$min_post_dawn < -30
  out[is_hour_early, "min_post_dawn"] = out[is_hour_early, "min_post_dawn"] + 60

  # Site 35031 is off by almost exactly an hour, & is in Indiana, which has
  # messy timezones/DST https://en.wikipedia.org/wiki/Time_in_Indiana
  is_hour_late = out$site_id == 35031 & out$min_post_dawn > 30
  out[is_hour_late, "min_post_dawn"] = out[is_hour_late, "min_post_dawn"] - 60

  out
}
