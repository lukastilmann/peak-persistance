# Function to count peaks in a base function
count_peaks_in_function <- function(func, t_grid) {
  # Convert function to tidyfun object
  f_tf <- tfd(func(t_grid), t_grid)

  # Calculate derivatives
  slope <- tf_derive(f_tf, order = 1)

  # Find peaks (local maxima) - where slope changes from positive to negative
  peaks <- tf_where(
    slope,
    value < 0 & dplyr::lag(value, 1, value[1]) > 0
  )

  # Return number of peaks (length of the first element of the list)
  return(length(peaks[[1]]))
}
