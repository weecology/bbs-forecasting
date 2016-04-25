library(ggplot2)
library(viridis)
library(maps)

# Assumes "results" is a data frame with long, lat, and p (which determines
# intensity)
# Shear amount is how much to move the top to the right (so negative numbers
# move to the left)
make_shear_map = function(results, shear_amount = -1){
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

  ggplot(results, aes(x = x, y = y, color = p)) +
    geom_polygon(aes_(~x, ~y), data = corners, inherit.aes = FALSE, color = "black", fill = "white") +
    shear_borders("state") +
    geom_point(size = 2) +
    scale_color_viridis(guide = FALSE) +
    coord_equal() +
    theme_void()
}

