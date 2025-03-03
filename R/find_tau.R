#' Calculate Significance Threshold for Peak Curvature in Function Curves
#'
#' @description
#' This function identifies a threshold value (tau) above which the curvature of peaks
#' in function curves are considered significant. The process involves smoothing the
#' input curves, calculating derivatives to identify peaks, and then determining the
#' normalized curvature values at these peaks.
#'
#' @param function_curves A functional data object containing the curves to analyze.
#' @param time_grid A numeric vector representing the time points or grid on which the functions are evaluated.
#' @param percentile A numeric value between 0 and 100 specifying the percentile to use for threshold determination.
#' @param smoothing_f A numeric value between 0 and 1 specifying the smoothing parameter for the lowess function.
#'        Default is 0.05 (lower values = less smoothing, higher values = more smoothing).
#'
#' @return A numeric value representing the threshold (tau) at the specified percentile of normalized curvature values.
#'
#' @details
#' The function performs the following steps:
#' 1. Smooths the input function curves using lowess smoothing
#' 2. Calculates first and second derivatives
#' 3. Identifies peaks as points where the slope changes from positive to negative
#' 4. Calculates normalized curvature values at these peaks
#' 5. Returns the value at the specified percentile of all curvature values
#'
#' Note that curvature values are negated and normalized by the maximum absolute value,
#' with negative values truncated to zero.
#'
#' @examples
#' \dontrun{
#' # Assuming functions_data is a functional data object with time_grid
#' tau <- find_tau(functions_data, time_grid, percentile = 95, smoothing_f = 0.05)
#' }
#'
#' @importFrom tf tf_smooth tf_derive tf_where tf_domain
#' @importFrom dplyr lag
#' @importFrom stats quantile
#'
#' @export
find_tau <- function(function_curves, time_grid, percentile, smoothing_f = 0.05) {
  # Input validation
  if (!is.numeric(percentile) || percentile < 0 || percentile > 100) {
    stop("'percentile' must be a numeric value between 0 and 100")
  }
  # TODO:  nicer: checkmate::assert_number(percentile, lower = 0, upper = 100)

  if (!is.numeric(smoothing_f) || smoothing_f <= 0 || smoothing_f >= 1) {
    stop("'smoothing_f' must be a numeric value between 0 and 1")
  }
  # TODO: nicer: checkmate::assert_number(smoothing_f, lower = 0, upper = 1)
  #   also: not really necessary, tf_smooth will check its inputs anyway...

  # Smoothing
  function_curves <- tf::tf_smooth(function_curves, "lowess", f = smoothing_f)
  # TODO:  default smoothing_f seems low?

  # Calculating first and second derivative
  slope <- tf::tf_derive(function_curves, order = 1)
  second_derivative <- tf::tf_derive(function_curves, order = 2)

  # Find peaks (local maxima) - where slope changes from positive to negative
  peak_locs <- tf::tf_where(
    slope,
    value < 0 & dplyr::lag(value, 1, value[1]) > 0)

  # Check if any peaks were found
  if (length(unlist(peak_locs)) == 0) {
    warning("No peaks detected in the function curves")
    return(0)  # Return a default value when no peaks are found
  }

  # Filter out arg values which are not in second derivative, as approximation
  # method shrinks domain
  domain_min <- tf::tf_domain(second_derivative)[1]
  domain_max <- tf::tf_domain(second_derivative)[2]
  peak_locs <- lapply(peak_locs, function(inner_list) {
    numeric_values <- as.numeric(unlist(inner_list))
    numeric_values[numeric_values >= domain_min & numeric_values <= domain_max]
  })

  # Remove empty lists after filtering
  peak_locs <- peak_locs[sapply(peak_locs, length) > 0]

  # Check if any peaks remain after domain filtering
  if (length(unlist(peak_locs)) == 0) {
    warning("No peaks remain after filtering for derivative domain")
    return(0)  # Return a default value when no valid peaks are found
  }

  # Calculating normalized curvature values
  peak_second_deriv <- second_derivative[, peak_locs, matrix = FALSE]
  curvature_values <- lapply(peak_second_deriv, function(mat) {
    # Extract the 'value' column (assuming it's the second column)
    values <- mat[, "value"]

    # Multiply by -1 and normalize by the maximum absolute value
    values_neg <- -1 * values
    values_neg <- ifelse(values_neg < 0, 0, values_neg)

    max_abs_value <- max(abs(values_neg))
    if (max_abs_value > 0) {
      normalized_values <- values_neg / max_abs_value
    } else {
      normalized_values <- values_neg  # Avoid division by zero
    }

    return(normalized_values)
  })

  # Combine all processed values into a single vector and get value at specified
  # percentile
  curv_vector <- unlist(curvature_values, use.names = FALSE)

  # Final check for valid curvature values
  if (length(curv_vector) == 0) {
    warning("No valid curvature values calculated")
    return(0)
  }

  tau_percentile <- quantile(curv_vector, percentile * 0.01)
  return(tau_percentile)
}
