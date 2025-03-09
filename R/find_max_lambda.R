#' Find First Lambda Value Where Frobenius Norm of Difference
#' Between Inelastic and Elastic Mean Drops Below Threshold
#'
#' @description
#' This function determines the optimal regularization parameter (lambda) for
#' function alignment by finding the lowest lambda value where the Frobenius norm
#' of the difference between inelastic and elastic mean functions drops below a
#' specified threshold. It employs a binary search (divide and conquer) approach
#' to efficiently identify this value.
#'
#' @param function_curves An object containing functional data. Can be either:
#'        (1) a data frame containing a 'curves' column that is a tidyfun 'tfd' object, or
#'        (2) a tidyfun 'tfd' object directly.
#' @param start_val Numeric. Initial lambda value to test (default = 2).
#' @param threshold Numeric. The threshold value for the Frobenius norm below which
#'        a lambda value is considered acceptable (default = 1e-2).
#' @param min_bound Numeric. Minimum distance between upper and lower bounds for
#'        search termination (default = 1e-3).
#' @param max_iter Integer. Maximum number of iterations for the alignment algorithm
#'        (default = 10). Will be doubled if initial convergence fails.
#' @param penalty character. Specifies penalty used in function alignment.
#' @param max_search_steps Integer. Maximum number of binary search steps to perform
#'        (default = 20).
#' @param parallel Logical. Whether to use parallel processing for alignment
#'        (default = FALSE).
#'
#' @return Numeric. The lowest lambda value where the Frobenius norm of the
#'         difference between the inelastic and elastic mean functions is below
#'         the specified threshold.
#'
#' @details
#' The function first checks if the provided starting lambda value yields a
#' Frobenius norm below the threshold. If not, it searches for an upper bound
#' by doubling lambda until finding a value that satisfies the threshold.
#'
#' Once appropriate bounds are established, a binary search is conducted to refine
#' the lambda value, seeking the lowest value that still satisfies the threshold
#' criterion. This approach efficiently navigates the trade-off between alignment
#' quality and regularization strength.
#'
#' The Frobenius norm is normalized by the square root of the sum of squared
#' deviations of the inelastic mean from its average value, providing a relative
#' measure of improvement.
#'
#' If the alignment algorithm fails to converge with the initial max_iter value,
#' the function will automatically double max_iter and retry once before warning
#' the user.
#'
#' @examples
#' \dontrun{
#' library(tidyfun)
#'
#' # Create synthetic functional data
#' t_grid <- seq(0, 1, length.out = 100)
#' curves_data <- replicate(10, sin(2 * pi * t_grid + rnorm(1, 0, 0.5)))
#' tf_curves <- tf_spline(curves_data, arg = t_grid)
#'
#' # Direct tfd object input
#' optimal_lambda1 <- find_max_lambda(
#'   function_curves = tf_curves,
#'   start_val = 2,
#'   threshold = 0.01
#' )
#'
#' # In a data frame
#' data_with_curves <- data.frame(id = 1)
#' data_with_curves$curves <- tf_curves
#' optimal_lambda2 <- find_max_lambda(
#'   function_curves = data_with_curves,
#'   start_val = 2,
#'   threshold = 0.01
#' )
#'
#' print(paste("Optimal lambda value:", optimal_lambda1))
#' }
#'
#' @export
find_max_lambda <- function(function_curves, start_val = 2, threshold = 1e-2,
                            min_bound = 1e-3, max_iter = 10,
                            penalty = "roughness",
                            max_search_steps = 20, parallel = FALSE) {
  # Input validation
  if (is.data.frame(function_curves)) {
    if (!("curves" %in% names(function_curves))) {
      stop("When 'function_curves' is a data frame, it must contain a 'curves' column")
    }
    if (!inherits(function_curves$curves, "tfd")) {
      stop("The 'curves' column in 'function_curves' must be a tidyfun 'tfd' object")
    }
    curves <- function_curves$curves
  } else if (inherits(function_curves, "tfd")) {
    curves <- function_curves
  } else {
    stop("'function_curves' must be either a data frame with a 'tfd' object in the 'curves' column or a 'tfd' object directly")
  }

  # Extract t_grid from the functional data
  t_grid <- tf::tf_arg(curves)

  # Function doubles max_iter if no convergence for initial value
  max_iter_flex <- max_iter

  # Calculate Frobenius norm for a given lambda
  calc_frob_norm <- function(lambda, norm_factor) {
    aligned_res <- align_functions(function_curves,
                                   lambda = lambda,
                                   max_iter = max_iter_flex,
                                   penalty = penalty,
                                   parallel = parallel,
                                   t_grid = t_grid)

    aligned <- aligned_res$aligned_curves

    if (!aligned_res$converged) {
      if (max_iter_flex == max_iter) {
        message("Doubling max_iter for lambda search and trying again")
        max_iter_flex <- max_iter * 2
        aligned_res <- align_functions(function_curves,
                                       lambda = lambda,
                                       max_iter = max_iter_flex,
                                       penalty = penalty,
                                       parallel = parallel,
                                       t_grid = t_grid)
        aligned <- aligned_res$aligned_curves
      } else {
        message(sprintf("No convergence after %s iterations. Set max_iter higher.",
                        max_iter_flex))
      }
    }

    sqrd_diff <- (mean_unaligned - mean(aligned))^2
    sqrt(sum(sqrd_diff[, t_grid])) * norm_factor
  }

  # Compute inelastic mean and Frobenius norm
  mean_unaligned <- mean(curves)
  mean_unaligned_vec <- mean_unaligned[, t_grid]
  squared_deviation <- (mean(mean_unaligned_vec) - mean_unaligned_vec)^2
  normalizing_factor <- 1 / sqrt(sum(squared_deviation))  # Invert to multiply later

  # Evaluate at start value
  fn_val_start <- calc_frob_norm(start_val, normalizing_factor)

  # If start value is not below threshold, search for upper bound
  if (fn_val_start > threshold) {
    start_i <- 1
    lower_bound <- 0
    fn_val_lower_bound <- 1
    candidate_fn <- fn_val_start
    candidate <- start_val

    while (candidate_fn > threshold) {
      # Get new candidate by extrapolating from prev lower bound and evaluated
      # candidate
      fn_lower_cand_diff <- fn_val_lower_bound - candidate_fn
      fn_lower_target_diff <- fn_val_lower_bound - threshold
      diff_ratio <- (fn_lower_target_diff / fn_lower_cand_diff) - 1
      interval_diff <- candidate - lower_bound

      # Setting lower bound to old candidate
      fn_val_lower_bound <- candidate_fn
      lower_bound <- candidate

      # Evaluating new candidate
      candidate_new <- candidate + (diff_ratio ^ (1/3)) * interval_diff
      candidate_fn <- calc_frob_norm(candidate_new, normalizing_factor)
      candidate <- candidate_new

      start_i <- start_i + 1

      if (start_i > max_search_steps / 2) {
        warning("Could not find lambda with norm below threshold. Consider increasing start_val or threshold.")
        return(upper_bound)  # Return the best found so far
      }
    }
    # Valid upper bound has been found
    upper_bound <- candidate_new
    fn_val_upper_bound <- candidate_fn
    lowest_valid <- candidate_new

  } else {
    # Start value was below threshold, can start divide and conquer approach
    lower_bound <- 0
    fn_val_lower_bound <- 1

    upper_bound <- start_val
    fn_val_upper_bound <- fn_val_start
    lowest_valid <- upper_bound

    start_i <- 1
  }

  # Divide and conquer search
  for (i in start_i:max_search_steps) {
    # Find new candidate by interpolating between bounds
    # Assumption that curve is roughly proportional to cubic root function
    fn_bound_diff <- fn_val_lower_bound - fn_val_upper_bound
    fn_target_lower_diff <- fn_val_lower_bound - threshold
    diff_ratio <- fn_target_lower_diff / fn_bound_diff
    bound_diff <- upper_bound - lower_bound

    # Making sure step as least as large as min bound
    if(diff_ratio^3 < 0.5){
      dist_from_lower <- max(diff_ratio^3 * bound_diff, min_bound)
      candidate_new <- lower_bound + dist_from_lower
    } else {
      dist_from_upper <- max((1 - diff_ratio^3) * bound_diff, min_bound)
      candidate_new <- upper_bound - dist_from_upper
    }

    # Evaluate candidate
    candidate_fn <- calc_frob_norm(candidate_new, normalizing_factor)

    # Check if evaluated value below threshold, set new lowest valid solution if so
    if (candidate_fn <= threshold){
      upper_bound <- candidate_new
      lowest_valid <- candidate_new
      fn_val_upper_bound <- candidate_fn
    } else {
      lower_bound <- candidate_new
      fn_val_lower_bound <- candidate_fn
    }

    # Stop if bounds are very close
    if (abs(lower_bound - upper_bound) / min_bound <= 1.001) {
      break
    }
  }

  return(lowest_valid)
}
