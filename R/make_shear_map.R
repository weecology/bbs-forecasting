#' @param results a data frame with long, lat, and p (which color determines 
#' intensity and is between 0 and 1). Note that any points with p==0 will not 
#' be plotted
#' @param shear_amout how much to move the top to the right (so negative numbers
#' move to the left)
#' @param map_color color for landmass in the background
#' @export
#' @importFrom ggplot2 geom_polygon geom_point
#' @importFrom ggplot2 ggplot aes aes_
#' @importFrom ggplot2 map_data
#' @importFrom ggplot2 theme_void coord_equal
#' @importFrom viridis scale_color_viridis
make_shear_map = function(results, shear_amount = -1, map_color = "#FEFFEF"){
  shear = matrix(c(1, shear_amount, 0, 1), nrow = 2)
  
  results[ , c("x", "y")] = as.matrix(results[,c("long", "lat")]) %*% shear
  
  # Sheared version of the US state map
  shear_borders = function (database = "world", regions = ".", fill = NA, colour = "grey50",
                            xlim = NULL, ylim = NULL, ...) {
    df <- map_data(database, regions, xlim = xlim, ylim = ylim)
    df[ , c("long", "lat")] = as.matrix(df[ , c("long", "lat")]) %*% shear
    geom_polygon(aes_(~long, ~lat, group = ~group), data = df,
                 fill = fill, colour = colour, ..., inherit.aes = FALSE)
  }
  
  # Define the parallelogram that will represent the map area
  xlim = range(map_data("state")$long)
  ylim = range(map_data("state")$lat)
  corner_matrix = cbind(x = rep(xlim, each = 2), y = c(ylim, rev(ylim)))
  corners = corners = structure(as.data.frame(corner_matrix %*% shear), names = c("x", "y"))
  
  results$alpha = ifelse(results$p == 0, 0, 1)
  
  ggplot(results, aes(x = x, y = y, color = p)) +
    geom_polygon(aes_(~x, ~y), data = corners, inherit.aes = FALSE, color = "black", fill = "white") +
    shear_borders("state", fill = map_color) +
    geom_point(size = 2, alpha = results$alpha) +
    scale_color_viridis(guide = FALSE, option = "magma", limits = c(0, 1),
                        direction = -1) +
    coord_equal() +
    theme_void()
}

