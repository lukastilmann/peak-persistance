#' Analyze the Mean Function of Functional Data
#'
#' @description
#' This function computes the mean of a collection of functional data objects and
#' analyzes its key features including peaks (local maxima), valleys (local minima),
#' peak heights, and peak curvatures. It uses derivatives to identify critical points
#' and characterize the shape of the mean function.
#'
#' @param curves A `tfd` or `tfb` object representing a collection of functional
#'   data curves. Must inherit from the `tf` class in the tidyfun package.
#'
#' @return A list containing:
#'   \item{peaks}{Numeric vector with the locations (x-coordinates) of peaks in the mean function}
#'   \item{valleys}{Numeric vector with the locations (x-coordinates) of valleys in the mean function}
#'   \item{function_mean}{The mean function as a tidyfun object}
#'   \item{peak_heights}{Numeric vector with the heights (y-coordinates) of peaks}
#'   \item{peak_curvatures}{Numeric vector with normalized curvatures at peaks (between 0 and 1)}
#'
#' @details
#' The function identifies peaks and valleys by examining where the first derivative
#' changes sign. Peaks occur when the slope changes from positive to negative, and valleys
#' occur when the slope changes from negative to positive. The curvature at peaks is
#' calculated using the negative of the second derivative and is normalized to have a
#' maximum value of 1. Negative curvatures (which shouldn't occur at true peaks) are
#' set to zero.
#'
#' @note
#' Peak locations near the boundaries of the function domain may be excluded when
#' calculating curvature due to limitations in derivative calculations at boundaries.
#'
#' @examples
#' \dontrun{
#' library(tidyfun)
#' # Create a sample tidyfun dataset
#' x <- seq(0, 1, length.out = 100)
#' curves_data <- replicate(10, sin(2 * pi * x + rnorm(1, 0, 0.2)))
#' curves <- tf_spline(curves_data, arg = x)
#'
#' # Analyze the mean function
#' results <- analyze_mean_function(curves)
#'
#' # Plot the mean function and mark peaks
#' plot(results$function_mean)
#' points(results$peaks, results$peak_heights, col = "red", pch = 16)
#' }
#'
#' @export
analyze_mean_function <- function(curves) {
  # Input validation
  if (!inherits(curves, "tf")) {
    stop("Input 'curves' must be a tidyfun 'tf' object. Please convert your data using appropriate tidyfun functions.")
  }

  # Calculate mean function
  curves_mean <- mean(curves)

  # Calculate derivatives
  slope <- tf::tf_derive(curves_mean, order = 1)
  second_derivative <- tf::tf_derive(curves_mean, order = 2)

  # Find peaks (local maxima) - where slope changes from positive to negative
  peaks <- tf::tf_where(
    slope,
    value < 0 & dplyr::lag(value, 1, value[1]) > 0)

  # Find valleys (local minima) - where slope changes from negative to positive
  valleys <- tf::tf_where(
    slope,
    value > 0 & dplyr::lag(value, 1, value[1]) < 0)

  # Compute height and negative curvature at peak locations
  # Accounting for reduction of argument domain at edges of interval when calculating derivative
  peak_loc <- peaks[[1]]

  # Filter peak locations to be within the valid domain of the second derivative
  peak_loc <- peak_loc[peak_loc < max(tf::tf_arg(second_derivative)) &
                         peak_loc > min(tf::tf_arg(second_derivative))]

  # Extract peak heights from the mean function
  peak_heights <- curves_mean[, peak_loc]

  # Extract and negate curvatures at peak locations (negative second derivative)
  peak_curvatures <- -second_derivative[, peak_loc]

  # Ensure non-negative curvatures (true peaks should have negative second derivatives)
  peak_curvatures <- ifelse(peak_curvatures < 0, 0, peak_curvatures)

  # Normalize curvatures to have a maximum of 1
  max_curvature <- max(peak_curvatures)
  if (max_curvature > 0) {
    peak_curvatures <- peak_curvatures / max_curvature
  }

  # Return results as a list
  return(list(
    peaks = peak_loc,
    valleys = valleys[[1]],
    function_mean = curves_mean,
    peak_heights = peak_heights,
    peak_curvatures = peak_curvatures
  ))
}
