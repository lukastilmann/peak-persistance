
#' Find Maximum Lambda Value
#'
#' @param function_curves Function values matrix/vector
#' @param t_grid Time vector
#' @param max_val Maximum value for lambda search
#' @param n_grid Number of grid points
#' @param parallel Logical indicating whether to use parallel processing
#' @return Maximum lambda value
find_max_lambda <- function(function_curves, t_grid, max_val, n_grid, parallel = FALSE) {
  # Create sequence of lambda values for testing
  lam_grid <- seq(0, max_val, length.out = n_grid)

  # Get aligned functions
  aligned <- lapply(lam_grid, function(lam){
    fal <- align_functions(function_curves, lambda = lam, parallel = parallel,
                           t_grid, func = "tw")
  })

  # Frobenius norm of difference between aligned mean and unaligned mean
  frob_norms <- lapply(aligned, function(fun_aligned){
    sqrd_diff <- (mean(function_curves$curves) - mean(fun_aligned))^2
    # tidyfun for this?
    sqrt(sum(sqrd_diff[, t_grid]))
  })

  # Find stopping index where norm exceeds threshold
  idx_last <- sum(frob_norms > 1e-2)

  # Return maximum lambda value
  max_lam <- lam_grid[idx_last]
  return(max_lam)
}
