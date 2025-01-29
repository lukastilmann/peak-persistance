
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
find_max_lambda <- function(function_curves, t_grid, max_val, threshold = 1e-2,
                            max_iter = 20, min_bound = 1e-3, parallel = FALSE) {

  mean_unaligned <- mean(function_curves$curves)
  # Calculate Frobenius norm for a given lambda
  calc_frob_norm <- function(lambda) {
    aligned <- align_functions(function_curves, lambda = lambda,
                               parallel = parallel, t_grid, func = "tw")
    sqrd_diff <- (mean_unaligned - mean(aligned))^2
    sqrt(sum(sqrd_diff[, t_grid]))
  }

  # Initialize binary search bounds
  left <- 0
  right <- max_val
  lowest_valid <- right  # Initialize to maximum value

  # Binary search
  for (i in 1:max_iter) {
    mid <- (left + right) / 2
    norm_val <- calc_frob_norm(mid)

    if (norm_val <= threshold) {
      # Found a valid lambda, but might not be the lowest
      # Store it and search lower
      right <- mid
      lowest_valid <- mid
    } else {
      # Current lambda is too small, search higher
      left <- mid
    }

    # Stop if bounds are very close
    if (abs(right - left) < min_bound) {
      break
    }
  }

  return(lowest_valid)
}
