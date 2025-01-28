library(tidyverse)
library(tidyfun)

# Computes mean function, locations of peaks and valleys, and height and
# curvature of peaks
analyze_mean_function <- function(curves) {

  curves_mean <- mean(curves)
  #TODO: vary smoothing by arguments to function
  #curves_mean <- tf_smooth(curves_mean)
  slope <- tf_derive(curves_mean, order = 1)
  second_derivative <- tf_derive(curves_mean, order = 2)

  # Find peaks (local maxima) - where slope changes from positive to negative
  peaks <- tf_where(
    slope,
    value < 0 & dplyr::lag(value, 1, value[1]) > 0)

  # Find valleys (local minima) - where slope changes from negative to positive
  valleys <- tf_where(
    slope,
    value > 0 & dplyr::lag(value, 1, value[1]) < 0)

  # Compute height and negative curvature at peak locations
  # Accounting for reduction of arg at edges of interval when calculating
  # derivative
  peak_loc <- peaks[[1]]
  peak_loc <- peak_loc[peak_loc < max(tf_arg(second_derivative)) &
                         peak_loc > min(tf_arg(second_derivative))]
  peak_heights <- curves_mean[, peak_loc]
  peak_curvatures <- - second_derivative[, peak_loc]

  # ensuring non-negative curvatures and normalizing
  # is there a way to do this in tidyfun?
  peak_curvatures <- ifelse(peak_curvatures < 0, 0, peak_curvatures)
  max_curvature <- max(peak_curvatures)

  if (max_curvature > 0){
    peak_curvatures <- peak_curvatures / max_curvature
  }

  return(list(
    peaks = peaks[[1]],
    valleys = valleys[[1]],
    function_mean = curves_mean,
    peak_heights = peak_heights,
    peak_curvatures = peak_curvatures
  ))
}
