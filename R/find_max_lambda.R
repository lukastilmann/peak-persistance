
#' Find First Lambda Value Where Frobenius Norm of Difference
#' Between Inelastic and Elastic Mean Drops Below Threshold
#'
#' Uses binary search to find the lowest lambda value where the Frobenius norm
#' first goes below the threshold
#'
#' @param function_curves List containing curves data
#' @param t_grid Time grid points
#' @param max_val Maximum value for lambda search
#' @param threshold Threshold for Frobenius norm (default: 1e-2)
#' @param min_bound Minimum for bound between upper and lower bound for search
#' @param max_iter Maximum number of binary search iterations (default: 20)
#' @param parallel Logical indicating whether to use parallel processing
#' @return Lowest lambda value that brings Frobenius norm below threshold
find_max_lambda <- function(function_curves, t_grid, start_val, threshold = 1e-2,
                            min_bound = 1e-3, max_iter = 10,
                            max_search_steps = 20, parallel = FALSE) {

  # Function doubles max_iter if no convergence for initial value
  max_iter_curr <- max_iter


  # Calculate Frobenius norm for a given lambda
  calc_frob_norm <- function(lambda, norm_factor) {
    aligned_res <- align_functions(function_curves,
                               lambda = lambda,
                               max_iter = max_iter_curr,
                               parallel = parallel,
                               t_grid = t_grid)
    aligned <- aligned_res$aligned_curves
    if (!aligned_res$converged){
      if (max_iter_curr == max_iter){
        message("Doubling max_iter for lambda search and trying again")
        max_iter_curr <- max_iter * 2
        aligned_res <- align_functions(function_curves,
                                       lambda = lambda,
                                       max_iter = max_iter_curr,
                                       parallel = parallel,
                                       t_grid = t_grid)
        aligned_res$aligned_curves
      } else {
        message(sprintf("No convergence after %s iterations. Set max_iter higher.",
                        max_iter_current))
      }
    }
    sqrd_diff <- (mean_unaligned - mean(aligned))^2
    sqrt(sum(sqrd_diff[, t_grid])) * norm_factor
  }

  # Compute inelastic mean and Frobenius norm of
  mean_unaligned <- mean(function_curves$curves)
  mean_unaligned_vec <- mean_unaligned[, t_grid]
  squared_dev <- (mean(mean_unaligned_vec) - mean_unaligned_vec)^2
  normalizing_factor <- sqrt(sum(squared_dev))

  # Evaluate at start value
  start_fn <- calc_frob_norm(start_val, normalizing_factor)
  if (start_fn > threshold){
    start_i = 1
    lower_bound <- start_val
    upper_fn <- start_fn
    upper_bound <- lower_bound
    # Search for upper bound
    while (upper_fn > threshold){
      upper_bound <- upper_bound * 2
      upper_fn <- calc_frob_norm(upper_bound, normalizing_factor)
      start_i = start_i + 1
    }
  } else {
    lower_bound <- 0
    upper_bound <- start_val
    lowest_valid <- upper_bound
    start_i = 1
  }

  # Divide and conquer search
  for (i in start_i:max_search_steps) {
    mid <- (upper_bound + lower_bound) / 2
    mid_fn <- calc_frob_norm(mid, normalizing_factor)

    if (mid_fn <= threshold) {
      # Found a valid lambda, but might not be the lowest
      # Store it and search lower
      upper_bound <- mid
      lowest_valid <- mid
    } else {
      # Current lambda is too small, search higher
      lower_bound <- mid
    }

    # Stop if bounds are very close
    if (abs(lower_bound - upper_bound) < min_bound) {
      break
    }
  }

  return(lowest_valid)
}
