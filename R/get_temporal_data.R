#' @importFrom rgdal readOGR
get_time_zones = function() {
  if(!file.exists("data/tz_world.zip")) {
    download.file("http://efele.net/maps/tz/world/tz_world.zip", "data/tz_world.zip")
    unzip("data/tz_world.zip", exdir = "data/tz_world")
  }
  readOGR("data/tz_world/world/" ,"tz_world", verbose = FALSE)
}

#' @importFrom rgeos gDistance
#' @importFrom sp proj4string SpatialPoints CRS over
add_time_zones = function(locations, plot_tz = FALSE){
  boundaries = get_time_zones()
  time_zone = over(
    SpatialPoints(locations[, c("long", "lat")],
                  proj4string = CRS(proj4string(boundaries))),
    boundaries
  )[[1]]

  out = cbind(locations, time_zone = time_zone)
  out$time_zone = as.character(out$time_zone)


  # locations of NAs
  missing_locations = out %>% filter(is.na(time_zone)) %>% distinct(long, lat, site_id) %>%
    SpatialPointsDataFrame(
      coords = .[1:2],
      data = .,
      coords.nrs = 1:2,
      proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    )

  missing_tzs = apply(
    gDistance(missing_locations, boundaries, byid=TRUE),
    2,
    function(x){as.character(boundaries@data[which.min(x), 1])}
  )

  for (i in 1:nrow(missing_locations)) {
    out[out$site_id == missing_locations$site_id[i], "time_zone"] = missing_tzs[i]
  }


  # Visualize the time zones to make sure that the imputed ones match their
  # neighbors
  if (plot_tz) {
    out %>%
      distinct(long, lat, time_zone, site_id) %>%
      mutate(missing = site_id %in% missing_locations$site_id) %>%
      plot(
        lat ~ long,
        data = .,
        col = factor(time_zone),
        cex = ifelse(missing, 2, .7),
        pch = as.integer(factor(time_zone)) %% 10,
        lwd = ifelse(missing, 3, .25))
  }


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
#' @importFrom rgdal readOGR
#' @importFrom lubridate day month yday year ymd_hm with_tz
#' @importFrom maptools crepuscule
get_temporal_data = function(start_yr, end_yr, min_num_yrs){

  events = get_bbs_data() %>%
    ungroup() %>%
    filter(site_id %in% get_bbs_gimms_ndvi()[["site_id"]]) %>%
    select(site_id, lat, long, year, month, day, start_time) %>%
    distinct() %>%
    add_time_zones() %>%
    group_by(time_zone) %>%
    do(add_times(.)) %>%
    ungroup() %>%
    mutate(yday = yday(date_time)) %>%
    dplyr::select(-time_zone, -month, -day)

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
