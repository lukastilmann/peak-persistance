
# Function to find template functional
template_function <- function(fn_mean, peak_locs, valley_locs, t_grid) {

  #TODO: only significant peaks and surrounding valleys
  # Beginning and end of grid, peaks and valleys
  idx <- sort(unique(c(t_grid[1], peak_locs, valley_locs, t_grid[length(t_grid)])))

  # Extract the corresponding y-values
  selected_points <- fn_mean[, idx]

  # Perform interpolation using splines (R's equivalent to MATLAB's pchip)
  ttemp <- splinefun(idx, selected_points, method = "monoH.FC")

  return(ttemp)
}
